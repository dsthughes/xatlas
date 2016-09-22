#include "Xatlas.hpp"
#include "Bam.hpp"
#include "CoverageCounter.hpp"
#include "EventHolder.hpp"
#include "IndelEvent.hpp"
#include "Sequence.hpp"
#include "SnpEvent.hpp"
#include "VcfWriter.hpp"
#include "faidx.h"
#include "sam.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <unistd.h>

/* ------===< xAtlas >===------ */

static const char *MY_NAME = "xatlas", *VERSION = "v0.0.2-hgsc";
static const char *OPTIONS = "\
Required arguments:\n\
    -r, --ref REF           Reference genome in FASTA format\n\
    -b, --bam BAM           Sorted and indexed input BAM file\n\
    -s, --sample-name SN    Sample name to use in the output VCF file\n\
    -p, --prefix PFX        Output VCF file prefix\n\
\n\
Options:\n\
    -c, --capture-bed BED       BED file of regions to process\n\
    -S, --no-scavenge           Equal majority SNP alt allele selection only by highest base quality\n\
    -a, --all-non-ref           Produce calls for all detected variants\n\
    -C, --no-cleanup            Do not remove temporary per-contig VCF files\n\
    -m, --min-snp-mapq MAPQ     Minimum read mapping quality for calling SNPs\n\
    -n, --min-indel-mapq MAPQ   Minimum read mapping quality for calling indels\n\
    -M, --max-coverage COV      Maximum coverage for calling variants normally\n\
    -x, --xindel                Use 6-parameter model trained on HiSeq X NA12878 data for indels\n\
    -g, --gvcf                  Include non-variant gVCF blocks in output VCF file\n\
    -A, --block-abs-lim LIM     gVCF non-variant block absolute range limit\n\
    -R, --block-rel-lim LIM     gVCF non-variant block relative range limit coefficient\n\
    -h, --help                  Show this help\n\
";

typedef struct args
{
    pthread_barrier_t barrier1;
    pthread_barrier_t barrier2;

    bam1_t **rec_buff;
    bam1_t **read_start;
    bam1_t **process_start;

    uint32_t read_section;
    uint32_t buff_size;
    uint32_t num_sections;
    uint32_t section_size;
    uint32_t process_size;

    samFile *sf;
    hts_itr_t *iter;

    EventHolder *events;
    Sequence *sequences;
    CoverageCounter *coverages;
    VcfWriter *writer;
    qual_t min_map_qual;
    qual_t min_map_qual2;

    bed_coord_list_t *segs;
    size_t idx;
    opts_s *opts;
    std::ofstream *oindelv;
    std::ofstream *osnpv;
} args_s;

void *read_bam(void *arg)
{
    args_s *args = (args_s *)arg;
    bool more_to_read = true;

    while (more_to_read)
    {
        uint32_t num_read = 0;
        bam1_t **curr_rec = args->read_start;

        while (num_read < args->section_size && sam_itr_next(args->sf, args->iter, *curr_rec) >= 0)
        {
            ++curr_rec;
            ++num_read;
        }

        if (num_read == 0)
            more_to_read = false;

        pthread_barrier_wait(&args->barrier1);

        args->process_start = args->read_start;
        args->process_size = num_read;

        args->read_section = (args->read_section + 1) % args->num_sections;
        if (args->read_section == 0)
            args->read_start = args->rec_buff;
        else
            args->read_start = curr_rec;

        pthread_barrier_wait(&args->barrier2);
    }

    //std::cerr << "[read_bam] exit\n";
    pthread_exit(nullptr);
}

void *process_records_indel(void *arg)
{
    const int16_t filtered_flags = BAM_FUNMAP + BAM_FSECONDARY + BAM_FQCFAIL + BAM_FDUP;
    args_s *args = (args_s *)arg;
    bool more_to_process = true, first_pass = true;
    uint32_t stagger = 50, num_processed = 0;
    bed_coord_t &seg = (*args->segs)[args->idx];

    while (more_to_process)
    {
        if (!first_pass)
        {
            bam1_t **next_rec = args->process_start;
            uint32_t num_read = 0;

            while (num_read < args->process_size)
            {
                bam1_t *curr_rec = *next_rec;
                ++next_rec;
                ++num_read;

                if ((curr_rec->core.flag & filtered_flags) != 0 ||
                    curr_rec->core.n_cigar == 0 ||
                    curr_rec->core.qual < args->min_map_qual ||
                    bam_aux2i(bam_aux_get(curr_rec, "NM")) == -1 ||
                    args->coverages->isToxic(curr_rec->core.pos, COVSIDX_INDEL_DP, curr_rec->core.l_qseq))
                {
                    continue;
                }

                args->events->collectIndels(curr_rec, *args->sequences, *args->coverages);
                ++num_processed;

                if (num_processed % stagger == 0)
                    args->writer->printIndelBuffer(curr_rec->core.pos, seg, *args->oindelv);
            }
        }

        pthread_barrier_wait(&args->barrier1);
        pthread_barrier_wait(&args->barrier2);

        first_pass = false;
        more_to_process = (args->process_size != 0);
    }

    // final variants and gvcf block
    args->writer->printIndelBuffer(args->opts->last_call_max, seg, *args->oindelv);

    if (args->opts->gvcftest)
    {
        if (args->events->last_call_indel < seg.first)
            args->events->last_call_indel = seg.first - 1;
        args->writer->outputGvcf(seg.second,
                                 args->events->last_call_indel,
                                 seg.second,
                                 *args->oindelv,
                                 QUALSIDX_INDEL_P,
                                 COVSIDX_INDEL_VR,
                                 COVSIDX_INDEL_RR,
                                 COVSIDX_INDEL_DP);
    }

    //std::cerr << "[process_reads_indel] exit\n";
    pthread_exit(nullptr);
}

void *process_records_snp(void *arg)
{
    const int16_t filtered_flags = BAM_FUNMAP + BAM_FSECONDARY + BAM_FQCFAIL + BAM_FDUP;
    args_s *args = (args_s *)arg;
    bool more_to_process = true, first_pass = true;
    uint32_t stagger = 10, num_processed = 0;
    bed_coord_t &seg = (*args->segs)[args->idx];

    while (more_to_process)
    {
        if (!first_pass)
        {
            bam1_t **next_rec = args->process_start;
            uint32_t num_read = 0;

            while (num_read < args->process_size)
            {
                bam1_t *curr_rec = *next_rec;
                ++next_rec;
                ++num_read;

                if ((curr_rec->core.flag & filtered_flags) != 0 ||
                    curr_rec->core.n_cigar == 0 ||
                    curr_rec->core.qual < args->min_map_qual ||
                    curr_rec->core.qual < args->min_map_qual2 ||
                    bam_aux2i(bam_aux_get(curr_rec, "NM")) == -1 ||
                    args->coverages->isToxic(curr_rec->core.pos, COVSIDX_SNP_DP, curr_rec->core.l_qseq))
                {
                    continue;
                }

                args->events->collectSnps(curr_rec, *args->sequences, *args->coverages);
                ++num_processed;

                if (num_processed % stagger == 0)
                    args->writer->printSnpBuffer(curr_rec->core.pos, seg, *args->osnpv);
            }
        }

        pthread_barrier_wait(&args->barrier1);
        pthread_barrier_wait(&args->barrier2);

        first_pass = false;
        more_to_process = (args->process_size != 0);
    }

    // final variants and gvcf block
    args->writer->printSnpBuffer(args->opts->last_call_max, seg, *args->osnpv);

    if (args->opts->gvcftest)
    {
        if (args->events->last_call_snp < seg.first)
            args->events->last_call_snp = seg.first - 1;
        args->writer->outputGvcf(seg.second,
                                 args->events->last_call_snp,
                                 seg.second,
                                 *args->osnpv,
                                 QUALSIDX_SNP_P,
                                 COVSIDX_SNP_VR,
                                 COVSIDX_SNP_RR,
                                 COVSIDX_SNP_DP);
    }

    //std::cerr << "[process_reads_snp] exit\n";
    pthread_exit(nullptr);
}

void do_a_thing(bed_coord_list_t &segs, const char *region, const char *bf, const char *refseq, const char *pfx, opts_s *opts, ill_s *ill, args_s *args)
{
    std::string basename(std::string(pfx) + "_contig_" + std::string(region));
    std::string indel_fn(basename + "_indel.vcf"), snp_fn(basename + "_snp.vcf");
    std::ofstream oindelv, osnpv;

    oindelv.open(indel_fn.c_str(), std::fstream::out | std::fstream::trunc);
    if (!oindelv)
    {
        std::cerr << "Unable to open output indel vcf file" << std::endl;
        exit(EXIT_FAILURE);
    }
    oindelv.precision(5);

    osnpv.open(snp_fn.c_str(), std::fstream::out | std::fstream::trunc);
    if (!osnpv)
    {
        std::cerr << "Unable to open output snp vcf file" << std::endl;
        exit(EXIT_FAILURE);
    }
    osnpv.precision(5);

    std::cerr << "Analysing contig " << region << std::endl;

    EventHolder events(opts->NEAR_END, ill->SLX_SNP_MAX_SUB, ill->SLX_SNP_MAX_INDEL);
    Sequence sequences(refseq, region);
    CoverageCounter coverages(sequences.len, ill->SLX_INDEL_LOW_MAP_QUAL_MAX_COVERAGE);
    VcfWriter writer(coverages, sequences, events, region, opts, ill);
    Bam bam(bf);
    std::vector< std::string > piece;

    std::cerr << "Fetched sequence for contig " << region << " of length " << sequences.len << std::endl;

    if (segs.empty() && opts->capturebed == nullptr)
    {
        piece.push_back(region);
        segs.push_back(std::make_pair(0, sequences.len - 1));
    }
    else
    {
        for (const auto &seg : segs)
        {
            if (seg.second >= seg.first)
            {
                std::stringstream x;
                x << region << ":" << (seg.first + 1) << "-" << (seg.second + 1);
                piece.push_back(x.str());
            }
        }
    }

    args->read_start = args->rec_buff;
    args->process_start = args->rec_buff;
    args->read_section = 0;
    args->process_size = 0;
    args->sf = bam.sf;
    args->events = &events;
    args->sequences = &sequences;
    args->coverages = &coverages;
    args->segs = &segs;
    args->writer = &writer;
    args->oindelv = &oindelv;
    args->osnpv = &osnpv;

    pthread_t read_bam_thread, process_records_snp_thread, process_records_indel_thread;
    void *status;

    for (size_t idx = 0; idx < piece.size(); ++idx)
    {
        std::cerr << "Processing reads for region " << region << ":" << (segs[idx].first + 1) << "-" << (segs[idx].second + 1) << std::endl;

        events.last_call_indel = segs[idx].first;
        events.last_call_snp = segs[idx].first;

        bam.setIter(piece[idx].c_str());
        args->iter = bam.iter;
        args->idx = idx;

        if (pthread_create(&read_bam_thread, nullptr, read_bam, args) != 0)
        {
            std::cerr << "create read_bam_thread error" << std::endl;
            exit(EXIT_FAILURE);
        }

        if (pthread_create(&process_records_indel_thread, nullptr, process_records_indel, args) != 0)
        {
            std::cerr << "create process_records_indel_thread bam thread error" << std::endl;
            exit(EXIT_FAILURE);
        }

        if (pthread_create(&process_records_snp_thread, nullptr, process_records_snp, args) != 0)
        {
            std::cerr << "create process_records_snp_thread bam thread error" << std::endl;
            exit(EXIT_FAILURE);
        }

        if (pthread_join(read_bam_thread, &status) != 0)
        {
            std::cerr << "join read_bam_thread error" << std::endl;
            exit(EXIT_FAILURE);
        }

        if (pthread_join(process_records_indel_thread, &status) != 0)
        {
            std::cerr << "join process_records_indel_thread thread error" << std::endl;
            exit(EXIT_FAILURE);
        }

        if (pthread_join(process_records_snp_thread, &status) != 0)
        {
            std::cerr << "join process_records_snp_thread thread error" << std::endl;
            exit(EXIT_FAILURE);
        }

        if (events.snps.size() + events.indels.size() > 0)
        {
            if (!events.snps.empty())
                std::cerr << "Error: " << events.snps.size() << " remaining unprocessed SNPs." << std::endl;
            if (!events.indels.empty())
                std::cerr << "Error: " << events.indels.size() << " remaining unprocessed INDELs." << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    oindelv.close();
    osnpv.close();
}

void do_all_the_things(bed_list_t bedlist, const char *bam, const char *refseq, const char *pfx, opts_s *opts, ill_s *ill, args_s *args)
{
    for (auto &seg : bedlist)
        do_a_thing(seg.second, seg.first.c_str(), bam, refseq, pfx, opts, ill, args);
}

std::string get_vcf_header(const char *sn, int argc, char **argv, bool snp, const char *ref, contigs_list_t &contigs, opts_s *opts, ill_s *ill)
{
    std::stringstream hs;
    time_t today = time(nullptr);
    tm now;
    localtime_r(&today, &now);
    char date[9];
    date[8] = '\0';
    sprintf(date, "%d%02d%02d", 1900 + now.tm_year, now.tm_mon + 1, now.tm_mday);
    hs << "##fileformat=VCFv4.1\n"
       << "##fileDate=" << date << "\n"
       << "##source=" << *argv << " " << VERSION << "\n"
       << "##reference=" << ref << "\n"
       << "##command=";
    for (int i = 0; i < argc; ++i)
    {
        hs << argv[i];
        if (i != argc - 1)
            hs << " ";
    }
    hs << "\n";

    // FOMRAT
    hs << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
       << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n"
       << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n"
       << "##FORMAT=<ID=RR,Number=1,Type=Integer,Description=\"Reference Read Depth\">\n"
       << "##FORMAT=<ID=VR,Number=1,Type=Integer,Description=\"Major Variant Read Depth\">\n";

    // FILTER
    if (snp)
    {
        hs << "##FILTER=<ID=low_snpqual,Description=\"SNP posterior probability is less than " << ill->SLX_SNP_CUTOFF << "\">\n"
           << "##FILTER=<ID=low_VariantReads,Description=\"Variant read depth is less than " << opts->SNP_MIN_COV << "\">\n"
           << "##FILTER=<ID=low_VariantRatio,Description=\"Variant read ratio is less than " << opts->SNP_HET_MIN << "\">\n"
           << "##FILTER=<ID=low_coverage,Description=\"Total coverage is less than " << ill->SLX_SNP_MIN_COVERAGE << "\">\n";
    }
    else
    {
        hs << "##FILTER=<ID=low_qual,Description=\"Indel posterior probability is less than " << opts->INDEL_PR_CUTOFF << "\">\n"
           << "##FILTER=<ID=low_VariantReads,Description=\"Variant read depth is less than " << opts->MIN_VAR_READS << "\">\n"
           << "##FILTER=<ID=low_VariantRatio,Description=\"Variant read ratio is less than " << opts->MIN_VAR_RATIO << "\">\n"
           << "##FILTER=<ID=low_coverage,Description=\"Total coverage is less than " << opts->MIN_DEPTH_COVERAGE << "\">\n";
    }
    if (snp)
        hs << "##FILTER=<ID=high_coverage,Description=\"Total coverage is greater than " << opts->SLX_SNP_HIGH_MAP_QUAL_MAX_COVERAGE << "\">\n"
           << "##FILTER=<ID=single_strand,Description=\"All variant reads are in a single strand direction\">\n";
    else
        hs << "##FILTER=<ID=read_end_ratio,Description=\"Ratio of variants reads within " << (opts->NEAR_END + 2) << "bp of read end is greater than " << opts->MAX_NEAR_READ_END_RATIO << "\">\n";
    hs << "##FILTER=<ID=No_data,Description=\"No valid reads on this site\">\n"
       << "##FILTER=<ID=No_var,Description=\"No valid variants reads on this site\">\n";

    // INFO
    if (snp)
        hs << "##INFO=<ID=P,Number=1,Type=Float,Description=\"SNP p-value\">\n";
    else
        hs << "##INFO=<ID=P,Number=1,Type=Float,Description=\"Indel p-value\">\n";
    if (snp)
        hs << "##INFO=<ID=equal_majority,Number=0,Type=Flag,Description=\"The called SNP has an equal number of reads indicating another variant call and base was chosen by highest summed quality score\">\n";
    if (opts->gvcftest)
    {
        hs << "##INFO=<ID=DPX,Number=4,Type=Float,Description=\"Minimum, maximum, average, and sample standard deviation of DP values within this non-variant block\">\n"
           << "##INFO=<ID=RRX,Number=4,Type=Float,Description=\"Minimum, maximum, average, and sample standard deviation of RR values within this non-variant block\">\n"
           << "##INFO=<ID=VRX,Number=4,Type=Float,Description=\"Minimum, maximum, average, and sample standard deviation of VR values within this non-variant block\">\n"
           << "##INFO=<ID=breakblock,Number=1,Type=String,Description=\"Reason why this non-variant block was terminated\">\n"
           << "##INFO=<ID=BLOCKAVG_" << opts->block_label << ",Number=0,Type=Flag,Description=\"Non-variant block compression scheme " << opts->block_label << "\">\n"
           << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of this non-variant block\">\n";
    }

    // contig
    for (const auto &contig : contigs)
        hs << "##contig=<ID=" << contig.first << ",length=" << contig.second << ">\n";

    hs << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sn << "\n";

    return hs.str();
}

void usage()
{
    std::cerr << MY_NAME << ", " << VERSION << ", Daniel S. T. Hughes, Jesse Farek" << std::endl
              << __FILE__ ", " __DATE__ ", " __TIME__ << std::endl
              << "g++-" << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__ << " (" << __cplusplus << ")" << std::endl
              << std::endl
              << "Complete reimplementation, merging, fixes, optimisation, extension, model" << std::endl
              << "update of AtlasSNP2 & AtlasINDEL2 procedures. See Challis et al., BMC" << std::endl
              << "Bioinformatics, 2012 13:8 & Shen et al., Genome Res. 2010. 20: 273-280. This" << std::endl
              << "Pre-Release of " << MY_NAME << " is available from https://github.com/dsthughes/xatlas.git" << std::endl
              << "and released under CREATIVE COMMONS, Attribution-NonCommercial-NoDerivs 3.0" << std::endl
              << "United States. See https://creativecommons.org/licenses/by-nc-nd/3.0/." << std::endl
              << std::endl
              << OPTIONS << std::endl;
}

int main(int argc, char **argv)
{
    char *ref = nullptr, *bam = nullptr, *pfx = nullptr, *sample_name = nullptr;
    qual_t min_snp_qual = 0, min_indel_qual = 1;
    int tmp_min_indel_qual = 0, tmp_min_snp_qual = 1;
    bool cleanup = true;

    // Options

    opts_s opts;

    opts.indeltest = false;
    opts.allnonref = false;
    opts.gvcftest = false;
    opts.scavenge = true;
    opts.NEAR_END = 3;
    opts.MIN_VAR_READS = 2;
    opts.STRAND_DIR_FILTER = false;
    opts.MIN_DEPTH_COVERAGE = 5;
    opts.MIN_VAR_RATIO = 0.06;
    opts.SNP_MIN_PR = 0.109;
    opts.INDEL_PR_CUTOFF = opts.SNP_MIN_PR;
    opts.MAX_NEAR_READ_END_RATIO = 0.8;
    opts.INDEL_DEPTH_CUTOFF = 5;
    opts.INDEL_HOM_VAR_CUTOFF = 0.6;
    opts.INDEL_HET_CUTOFF = 0.06;
    opts.SNP_HET_MIN = 0.1;
    opts.SNP_HET_MAX = 0.9;
    opts.SNP_STRAND_RATIO_CUTOFF = 0.01;
    opts.SNP_STRAND_TEST_COV_CUTOFF = 16;
    opts.SNP_MIN_COV = 2;
    opts.block_abs_lim = 3; // 10;
    opts.block_rel_lim = 0.3;
    opts.block_rel_min = 1.0;
    opts.last_call_max = INT32_MAX;
    opts.capturebed = nullptr;
    opts.SLX_SNP_HIGH_MAP_QUAL_MAX_COVERAGE = 8000;

    // Illumina model

    ill_s ill;

    ill.SLX_SNP_intercept = -9.088;
    ill.SLX_SNP_quality_score = 0.162;
    ill.SLX_SNP_NQS_pass = 1.645;
    ill.SLX_SNP_swap = 0.0;
    ill.SLX_SNP_rel_pos = 2.349;
    ill.SLX_SNP_MAX_SUB = 0.05;
    ill.SLX_SNP_MAX_INDEL = 0.05;
    ill.SLX_INDEL_LOW_MAP_QUAL_MAX_COVERAGE = 16383;
    ill.SLX_SNP_MIN_COVERAGE = 6;
    ill.prior_err_c = 0.1;
    ill.prior_snp_c = 0.9;
    ill.SLX_SNP_CUTOFF = 0.95;

    ill.SLX_err_prior_arr[0] = 0.991;
    ill.SLX_err_prior_arr[1] = 0.004;
    ill.SLX_err_prior_arr[2] = 0.002;
    ill.SLX_err_prior_arr[3] = 0.001;
    ill.SLX_err_prior_arr[4] = 0.001;
    ill.SLX_err_prior_arr[5] = 0.001;
    ill.SLX_err_prior_arr[6] = 0.000002;  // 0.000005;
    ill.SLX_err_prior_arr[7] = 0.000005;  // 0.000002;
    ill.SLX_err_prior_arr[8] = 0.0000001; // ?
    ill.SLX_err_prior_arr[9] = 0.000001;  // ?

    ill.SLX_snp_prior_arr[0] = 0.0000001;
    ill.SLX_snp_prior_arr[1] = 0.0000001;
    ill.SLX_snp_prior_arr[2] = 0.0000001;
    ill.SLX_snp_prior_arr[3] = 0.0000001;
    ill.SLX_snp_prior_arr[4] = 0.0000001;
    ill.SLX_snp_prior_arr[5] = 0.0000001;
    ill.SLX_snp_prior_arr[6] = 0.0000001;
    ill.SLX_snp_prior_arr[7] = 0.0000001;
    ill.SLX_snp_prior_arr[8] = 0.027;
    ill.SLX_snp_prior_arr[9] = 0.973;

    ill.SLX_INDEL_intercept = -20.5000;
    ill.SLX_INDEL_local_entropy = 3.39065;
    ill.SLX_INDEL_strand_dir = 3.02573;
    ill.SLX_INDEL_norm_var_sq = 0.32695;
    ill.SLX_INDEL_mean_avg_nqs = 0.37184;

    ill.WGS_INDEL_intercept = -14.95535266;
    ill.WGS_INDEL_simple_local_entropy = 0.16228492;
    ill.WGS_INDEL_strand_dir = 4.84302146;
    ill.WGS_INDEL_norm_var_square = 0.02302639;
    ill.WGS_INDEL_mean_avg_nqs = 0.18140485;
    ill.WGS_INDEL_mean_map_qual = 0.10577339;
    ill.WGS_INDEL_mean_var_rate = -62.03643929;

    // Runtime options
    // TODO options for calling only snps/only indels

    int c, optidx;
    const char *short_options = "rbspcSaCmnMARxgh";
    const char *all_options = "0r:b:s:p:c:SaCm:n:M:A:R:xgh";
    static struct option long_options[] = {
        {"ref", 1, nullptr, 0},            // r
        {"bam", 1, nullptr, 0},            // b
        {"sample-name", 1, nullptr, 0},    // s
        {"prefix", 1, nullptr, 0},         // p
        {"capture-bed", 1, nullptr, 0},    // c
        {"no-scavenge", 0, nullptr, 0},    // S
        {"all-non-ref", 0, nullptr, 0},    // a
        {"no-cleanup", 0, nullptr, 0},     // C
        {"min-snp-mapq", 1, nullptr, 0},   // m
        {"min-indel-mapq", 1, nullptr, 0}, // n
        {"max-coverage", 1, nullptr, 0},   // M
        {"block-abs-lim", 1, nullptr, 0},  // A
        {"block-rel-lim", 1, nullptr, 0},  // R
        {"xindel", 0, nullptr, 0},         // x
        {"gvcf", 0, nullptr, 0},           // g
        {"help", 0, nullptr, 0},           // h
        {nullptr, 0, nullptr, 0}};

    if (argc == 1)
    {
        usage();
        return EXIT_SUCCESS;
    }

    while ((c = getopt_long(argc, argv, all_options, long_options, &optidx)) != -1)
    {
        if (c == 0)
            c = short_options[optidx];

        switch (c)
        {
        case 'r': // ref
            ref = optarg;
            break;
        case 'b': // bam
            bam = optarg;
            break;
        case 's': // sample-name
            sample_name = optarg;
            break;
        case 'p': // prefix
            pfx = optarg;
            break;
        case 'c': // capture-bed
            opts.capturebed = optarg;
            break;
        case 'S': // no-scavenge
            opts.scavenge = false;
            break;
        case 'a': // all-non-ref
            opts.allnonref = true;
            break;
        case 'C': // no-cleanup
            cleanup = false;
            break;
        case 'm': // min-snp-mapq
            tmp_min_snp_qual = (qual_t)atoi(optarg);
            break;
        case 'n': // min-indel-mapq
            tmp_min_indel_qual = (qual_t)atoi(optarg);
            break;
        case 'M': // max-coverage
            opts.SLX_SNP_HIGH_MAP_QUAL_MAX_COVERAGE = atoi(optarg);
            break;
        case 'A': // block-abs-lim
            opts.block_abs_lim = atoi(optarg);
            break;
        case 'R': // block-rel-lim
            opts.block_rel_lim = atof(optarg);
            break;
        case 'x': // xindel
            opts.indeltest = true;
            break;
        case 'g': // gvcf
            opts.gvcftest = true;
            break;
        case 'h': // help
            usage();
            return EXIT_SUCCESS;
        default:
            exit(EXIT_FAILURE);
        }
    }

    snprintf(opts.block_label, 24, "min%dp%da", (int)(100 * opts.block_rel_lim), (int)opts.block_abs_lim);

    // Begin

    std::cerr << "Running " << *argv << std::endl;

    LOGIT(ill.SLX_INDEL_LOW_MAP_QUAL_MAX_COVERAGE);
    LOGIT(ill.SLX_SNP_MIN_COVERAGE);
    LOGIT(opts.block_abs_lim);
    LOGIT(opts.block_rel_lim);
    LOGIT(opts.block_label);
    LOGIT(opts.SLX_SNP_HIGH_MAP_QUAL_MAX_COVERAGE);

    if (opts.SLX_SNP_HIGH_MAP_QUAL_MAX_COVERAGE > 8191)
    {
        std::cerr << "Maximum coverage is too high" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (bam != nullptr)
        std::cerr << "Using BAM '" << bam << "'" << std::endl;
    else
    {
        std::cerr << "No input BAM given" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (ref != nullptr)
        std::cerr << "Using referebce '" << ref << "'" << std::endl;
    else
    {
        std::cerr << "No reference given" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (sample_name != nullptr)
        std::cerr << "Using sample name '" << sample_name << "'" << std::endl;
    else
    {
        std::cerr << "No sample name given" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (tmp_min_indel_qual < 0 || tmp_min_indel_qual > 60)
    {
        std::cerr << "Minimum Indel mapping quality must be in the range [0, 60]" << std::endl;
        exit(EXIT_FAILURE);
    }
    else if (tmp_min_snp_qual < 0 || tmp_min_snp_qual > 60)
    {
        std::cerr << "Minimum SNP mapping quality must be in the range [0, 60]" << std::endl;
        exit(EXIT_FAILURE);
    }
    else if (tmp_min_snp_qual < tmp_min_indel_qual)
    {
        std::cerr << "Setting lower minimum mapping qual for SNPs has not been implemented" << std::endl;
        exit(EXIT_FAILURE);
    }

    min_snp_qual = tmp_min_indel_qual;
    min_indel_qual = tmp_min_snp_qual;

    // Files

    if (pfx == nullptr)
    {
        std::cerr << "No filename prefix given" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::cerr << "Prefix: " << pfx << std::endl;
    std::cerr << "Minimum snp map qual: " << (short)min_snp_qual << std::endl;
    std::cerr << "Minimum indel map qual: " << (short)min_indel_qual << std::endl;
    std::cerr << "Using BAM header for input sequence list" << std::endl;

    samFile *in = sam_open_format(bam, "rb", 0);
    if (in == nullptr)
    {
        std::cerr << "Failed to load BAM file \"" << bam << "\"" << std::endl;
        exit(EXIT_FAILURE);
    }

    bam_hdr_t *header = sam_hdr_read(in);
    if (header == nullptr)
    {
        std::cerr << "Failed to load header for BAM file \"" << bam << "\"" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::cerr << "Found " << header->n_targets << " contigs in BAM" << std::endl;

    // Which regions

    contigs_list_t contigs;

    for (int i = 0; i < header->n_targets; ++i)
    {
        std::cerr << "Contig " << i << " : " << header->target_name[i] << " : " << header->target_len[i] << std::endl;
        contigs.push_back(std::make_pair(std::string(header->target_name[i]), header->target_len[i]));
    }

    bam_hdr_destroy(header);
    sam_close(in);

    bed_list_t bedlist;

    if (opts.capturebed != nullptr)
    {
        auto *bedmap = new bed_coord_map_t;
        auto *tmp = new std::string;

        std::cerr << "Using " << opts.capturebed << " for input sequence list" << std::endl;
        std::ifstream bf(opts.capturebed);
        if (!bf)
        {
            std::cerr << "Problem reading " << opts.capturebed << std::endl;
            exit(EXIT_FAILURE);
        }

        char tname[64];
        pos_t start, end;
        while (getline(bf, *tmp) != nullptr)
        {
            if ((sscanf(tmp->c_str(), "%63s %d %d", tname, &start, &end)) < 3)
                continue;
            (*bedmap)[std::string(tname)].push_back(std::make_pair(start - 1, end - 1));
        }

        bf.close();

        for (const auto &contig : contigs)
            if (bedmap->count(contig.first) > 0)
                bedlist.push_back(std::make_pair(contig.first, (*bedmap)[contig.first]));

        delete bedmap;
        delete tmp;
    }
    else
        for (const auto &contig : contigs)
            bedlist.push_back(std::make_pair(contig.first, std::vector< bed_coord_t >()));

    if (bedlist.empty())
    {
        std::cerr << "Nothing to do" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::cerr << "Will analyse " << bedlist.size() << (bedlist.size() > 1 ? " contigs" : " contig") << std::endl;

    // Setup

    args_s args;

    args.opts = &opts;
    args.num_sections = 5;
    args.section_size = 10000;
    args.min_map_qual = min_snp_qual;
    args.min_map_qual2 = min_indel_qual;
    args.buff_size = args.num_sections * args.section_size;
    args.rec_buff = new bam1_t *[args.buff_size];
    for (uint32_t k = 0; k < args.buff_size; ++k)
        args.rec_buff[k] = bam_init1();

    pthread_barrier_init(&args.barrier1, nullptr, 3);
    pthread_barrier_init(&args.barrier2, nullptr, 3);

    do_all_the_things(bedlist, bam, ref, pfx, &opts, &ill, &args);

    // Cleanup

    pthread_barrier_destroy(&args.barrier1);
    pthread_barrier_destroy(&args.barrier2);

    for (uint32_t k = 0; k < args.buff_size; ++k)
        bam_destroy1(args.rec_buff[k]);
    delete[] args.rec_buff;

    // Final VCFs
    // TODO split this into 2 threads for indel/snp final vcfs
    // TODO or write directly to final vcfs in process threads

    std::cerr << "Writing final VCF" << std::endl;

    std::string snp_header, indel_header, snp_fn, indel_fn;
    std::string pfx_str(pfx);
    std::string snp_final_fn(pfx_str + "_snp.vcf");
    std::string indel_final_fn(pfx_str + "_indel.vcf");

    // snp header
    FILE *snp_final_fp = fopen(snp_final_fn.c_str(), "w");
    if (snp_final_fp == nullptr)
    {
        perror(nullptr);
        exit(EXIT_FAILURE);
    }
    snp_header = get_vcf_header(sample_name, argc, argv, true, ref, contigs, &opts, &ill);
    fputs(snp_header.c_str(), snp_final_fp);

    // indel header
    FILE *indel_final_fp = fopen(indel_final_fn.c_str(), "w");
    if (indel_final_fp == nullptr)
    {
        perror(nullptr);
        exit(EXIT_FAILURE);
    }
    indel_header = get_vcf_header(sample_name, argc, argv, false, ref, contigs, &opts, &ill);
    fputs(indel_header.c_str(), indel_final_fp);

    auto *buffer = new char[BUFSIZ];

    // for each contig
    for (const auto &contig : bedlist)
    {
        // concatenate snp vcf
        snp_fn = pfx_str + "_contig_" + contig.first + "_snp.vcf";

        FILE *snp_fp = fopen(snp_fn.c_str(), "r");
        if (snp_fp == nullptr)
        {
            perror(nullptr);
            exit(EXIT_FAILURE);
        }

        while (fgets(buffer, BUFSIZ, snp_fp) != nullptr)
            fputs(buffer, snp_final_fp);

        fclose(snp_fp);

        // concatenate indel vcf
        indel_fn = pfx_str + "_contig_" + contig.first + "_indel.vcf";
        FILE *indel_fp = fopen(indel_fn.c_str(), "r");
        if (indel_fp == nullptr)
        {
            perror(nullptr);
            exit(EXIT_FAILURE);
        }

        while (fgets(buffer, BUFSIZ, indel_fp) != nullptr)
            fputs(buffer, indel_final_fp);

        fclose(indel_fp);
    }

    delete[] buffer;
    fclose(snp_final_fp);
    fclose(indel_final_fp);

    // Delete temp files

    if (cleanup)
    {
        std::string fn;
        for (const auto &fit : bedlist)
        {
            fn = std::string(pfx) + "_contig_" + fit.first + "_snp.vcf";
            if (remove(fn.c_str()) != 0)
                std::cerr << "Unable to remove \"" << fn << "\"" << std::endl;

            fn = std::string(pfx) + "_contig_" + fit.first + "_indel.vcf";
            if (remove(fn.c_str()) != 0)
                std::cerr << "Unable to remove \"" << fn << "\"" << std::endl;
        }
    }

    return EXIT_SUCCESS;
}
