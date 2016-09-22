#include "VcfWriter.hpp"
#include "BlockStats.hpp"
#include "Logit.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <map>
#include <sstream>
#include <vector>

VcfWriter::VcfWriter()
{
}

VcfWriter::VcfWriter(CoverageCounter &_coverages,
                     Sequence &_sequences,
                     EventHolder &_events,
                     const char *_region,
                     const opts_s *_opts,
                     const ill_s *_ill)
    : coverages(&_coverages),
      sequences(&_sequences),
      events(&_events),
      region(_region),
      opts(_opts),
      ill(_ill)
{
}

VcfWriter::VcfWriter(const VcfWriter &) = default;

VcfWriter::~VcfWriter() = default;

void VcfWriter::outputGvcf(pos_t next_var_pos,
                           pos_t last_call_pos,
                           pos_t region_end_pos,
                           std::ofstream &output,
                           quals_idx_e i_quals_,
                           covs_idx_e i_vr_,
                           covs_idx_e i_rr_,
                           covs_idx_e i_dp_)
{
    if (last_call_pos == this->opts->last_call_max)
        return;

    BlockStats block_q(this->opts->block_abs_lim, this->opts->block_rel_lim);
    BlockStats block_vr(this->opts->block_abs_lim, this->opts->block_rel_lim);
    BlockStats block_rr(this->opts->block_abs_lim, this->opts->block_rel_lim);
    BlockStats block_dp(this->opts->block_abs_lim, this->opts->block_rel_lim);
    pos_t prev_block = last_call_pos;
    bool endofregion = false;
    std::string reason;

    for (pos_t i = last_call_pos; i <= next_var_pos; ++i)
    {
        bool breakit = true;
        coverage_t vr_cov = this->coverages->getCoverage(i, i_vr_);
        coverage_t rr_cov = this->coverages->getCoverage(i, i_rr_);
        coverage_t dp_cov = this->coverages->getCoverage(i, i_dp_);

        if (i == region_end_pos)
        {
            reason = "endofregion";
            endofregion = true;
        }
        else if (i == next_var_pos)
            reason = "blocked";
        else if (dp_cov == 0 && block_dp.min > 0)
            reason = "cov2nocov";
        else if (block_dp.max == 0 && dp_cov > 0 && block_dp.k > 0)
            reason = "nocov2cov";
        else if (block_dp.breakBlock((double)dp_cov))
            reason = "dp";
        else if (block_rr.breakBlock((double)rr_cov))
            reason = "rr";
        else if (block_vr.breakBlock((double)vr_cov))
            reason = "vr";
        else
            breakit = false;

        if (breakit)
        {
            output << this->region << "\t"
                   << prev_block + 1 << "\t"
                   << ".\t"
                   << this->sequences->seq[prev_block] << "\t"
                   << ".\t"
                   << ".\t"
                   << "PASS\t"
                   << "END=" << (endofregion ? i + 1 : i) << ";"
                   << "BLOCKAVG_" << this->opts->block_label << ";"
                   << "PX="
                   << block_q.min << ","
                   << block_q.max << ","
                   << block_q.getMean() << ","
                   << block_q.getStdDev() << ";"
                   << "VRX="
                   << block_vr.min << ","
                   << block_vr.max << ","
                   << block_vr.getMean() << ","
                   << block_vr.getStdDev() << ";"
                   << "RRX="
                   << block_rr.min << ","
                   << block_rr.max << ","
                   << block_rr.getMean() << ","
                   << block_rr.getStdDev() << ";"
                   << "DPX="
                   << block_dp.min << ","
                   << block_dp.max << ","
                   << block_dp.getMean() << ","
                   << block_dp.getStdDev() << ";"
                   << "breakblock=" << reason << "\t"
                   << "GT:VR:RR:DP:GQ\t"
                   << "0/0:"
                   << block_vr.min << ":"
                   << block_rr.min << ":"
                   << block_dp.min << ":"
                   << ".\n";

            if (endofregion)
                return;

            prev_block = i;
            block_q.resetBlocker();
            block_vr.resetBlocker();
            block_rr.resetBlocker();
            block_dp.resetBlocker();
        }

        block_q.addValue(this->coverages->getP(i, i_quals_));
        block_vr.addValue(vr_cov);
        block_rr.addValue(rr_cov);
        block_dp.addValue(dp_cov);
    }
}

inline void VcfWriter::addField(std::string &str, const char *to_add)
{
    if (!str.empty())
        str.push_back(';');
    str.append(to_add);
}

void VcfWriter::printSnpBuffer(pos_t next_var_pos, const bed_coord_t &seg, std::ofstream &osnpv)
{
    // TODO Remove binning step and fix cutoffs

    char genotype[4];
    std::string filter;
    std::map< char, uint32_t > counts;
    std::map< char, big_qual_t > quals;
    std::map< char, double > prs;
    SnpMap::iterator sv_it;

    for (sv_it = this->events->snps.begin(); sv_it != this->events->snps.end() && sv_it->first < next_var_pos; ++sv_it)
    {
        coverage_t refbase_cov = this->coverages->getCoverage(sv_it->first, COVSIDX_SNP_RR);
        coverage_t alternative_reads = sv_it->second.size();
        coverage_t total_coverage = this->coverages->getCoverage(sv_it->first, COVSIDX_SNP_DP) + alternative_reads;

        if (this->opts->allnonref || total_coverage > this->opts->SLX_SNP_HIGH_MAP_QUAL_MAX_COVERAGE)
            continue;

        char high_base = '.';
        coverage_t high_base_cov = 0;
        counts.clear();
        quals.clear();

        for (const auto &snp_it : sv_it->second)
        {
            ++counts[snp_it.allele_base];
            quals[snp_it.allele_base] += snp_it.qual;
        }

        for (const auto &counts_it : counts)
        {
            if (counts_it.second > high_base_cov)
            {
                high_base = counts_it.first;
                high_base_cov = counts_it.second;
            }
        }

        coverage_t equal_majority = 0;
        uint8_t bin;
        if (high_base_cov > 1)
        {
            for (const auto &counts_it : counts)
                if (counts_it.second == high_base_cov)
                    ++equal_majority;

            if (equal_majority >= 2)
            {
                if (!this->opts->scavenge)
                {
                    big_qual_t high_qual = quals[high_base];

                    for (const auto &counts_it : counts)
                        if (counts_it.second == high_base_cov && quals[counts_it.first] > high_qual)
                            high_base = counts_it.first;
                }
                else
                {
                    prs.clear();

                    for (const auto &snp_it : sv_it->second)
                    {
                        double pr = 1.0 - snp_logit(this->ill, (double)snp_it.qual, snp_it.nqs, (double)snp_it.type, snp_it.rel_pos);
                        if (prs.count(snp_it.allele_base) == 0)
                            prs[snp_it.allele_base] = pr;
                        else
                            prs[snp_it.allele_base] *= pr;
                    }

                    for (auto &prs_it : prs)
                    {
                        bin = floor((1.0 - prs_it.second) * 10);
                        if (bin > 9)
                            bin = 9;
                        prs_it.second = (1.0 / (1.0 + ((this->ill->SLX_err_prior_arr[bin] * this->ill->prior_err_c) /
                                                       (this->ill->SLX_snp_prior_arr[bin] * this->ill->prior_snp_c)))); ///
                    }

                    double high_pr = prs[high_base];

                    for (const auto &counts_it : counts)
                    {
                        if (counts_it.second == high_base_cov && prs[counts_it.first] > high_pr)
                        {
                            high_base = counts_it.first;
                            high_pr = prs[high_base]; ///
                        }
                    }
                }
            }
        }

        coverage_t var_cov = 0, pos_strand = 0;
        double pr_err_j_prod_sum = 1.0;

        for (const auto &snp_it : sv_it->second)
        {
            if (snp_it.allele_base == high_base)
            {
                ++var_cov;
                if (snp_it.read_strand)
                    ++pos_strand;
                /*
                pr_snp_i_read = snp_logit(this->ill, snp_it.qual, snp_it.nqs, snp_it.type, (double)snp_it.dist3 / snp_it.read_len);
                pr_err_i_read = 1 - pr_snp_i_read;
                */
                pr_err_j_prod_sum *= (1.0 - snp_logit(this->ill, (double)snp_it.qual, snp_it.nqs, (double)snp_it.type, snp_it.rel_pos));
            }
        }

        this->coverages->setCoverage(sv_it->first, COVSIDX_SNP_VR, var_cov);
        this->coverages->addCoverage(sv_it->first, COVSIDX_SNP_DP, var_cov);
        //this->coverages->set_force_call(COVSIDX_SNP_AR, sv_it->first, alternative_reads);

        //pr_snp_j_prod_sum = 1.0 - pr_err_j_prod_sum;
        bin = floor((1.0 - pr_err_j_prod_sum) * 10);
        if (bin > 9)
            bin = 9;

        /*
        pr_S_j_ERR_c = this->ill->SLX_err_prior_arr[bin];
        pr_S_j_SNP_c = this->ill->SLX_snp_prior_arr[bin];
        posterior_err = pr_S_j_ERR_c * this->ill->prior_err_c;
        posterior_snp = pr_S_j_SNP_c * this->ill->prior_snp_c;
        */
        double pr_SNP_S_j_c_j = 1.0 / (1.0 + ((this->ill->SLX_err_prior_arr[bin] * this->ill->prior_err_c) /
                                              (this->ill->SLX_snp_prior_arr[bin] * this->ill->prior_snp_c)));

        if (!this->opts->allnonref && pr_SNP_S_j_c_j < this->opts->SNP_MIN_PR)
            continue;

        if (this->opts->capturebed == nullptr ||
            (sv_it->first >= seg.first &&
             sv_it->first <= seg.second))
        {
            char refbase = this->sequences->seq[sv_it->first];
            coverage_t dp_cov = var_cov + refbase_cov;
            qual_t snpq = round(-10 * log10(1 - pr_SNP_S_j_c_j + 0.000001));
            filter.clear();

            if (pr_SNP_S_j_c_j < this->ill->SLX_SNP_CUTOFF)
                addField(filter, "low_snpqual");

            if (dp_cov == 0)
            {
                addField(filter, "No_data");
                strncpy(genotype, "./.", 4);
            }
            else
            {
                if (dp_cov < this->ill->SLX_SNP_MIN_COVERAGE)
                    addField(filter, "low_coverage");
                else if (dp_cov > this->opts->SLX_SNP_HIGH_MAP_QUAL_MAX_COVERAGE)
                    addField(filter, "high_coverage");

                double cov_ratio_pos = (double)pos_strand / var_cov;
                if (dp_cov >= this->opts->SNP_STRAND_TEST_COV_CUTOFF &&
                    std::min(cov_ratio_pos, 1.0 - cov_ratio_pos) < this->opts->SNP_STRAND_RATIO_CUTOFF)
                {
                    addField(filter, "single_strand");
                }

                if (var_cov == 0)
                    filter = "No_var";
                else if (var_cov <= this->opts->SNP_MIN_COV)
                    addField(filter, "low_VariantReads");

                double cov_ratio = (double)var_cov / dp_cov;
                if (cov_ratio <= this->opts->SNP_HET_MIN)
                {
                    strncpy(genotype, "0/0", 4);
                    addField(filter, "low_VariantRatio");
                }
                else if (cov_ratio < this->opts->SNP_HET_MAX)
                    strncpy(genotype, "0/1", 4);
                else
                    strncpy(genotype, "1/1", 4);
            }

            if (filter.empty())
                filter = "PASS";

            if (this->opts->gvcftest && this->events->last_call_snp < sv_it->first)
            {
                outputGvcf(sv_it->first,
                           this->events->last_call_snp,
                           seg.second,
                           osnpv,
                           QUALSIDX_SNP_P,
                           COVSIDX_SNP_VR,
                           COVSIDX_SNP_RR,
                           COVSIDX_SNP_DP);
            }
            this->events->last_call_snp = sv_it->first + 1;

            std::stringstream info_ss;
            info_ss << "P=" << pr_SNP_S_j_c_j;
            if (equal_majority >= 2)
                info_ss << ";equal_majority";

            osnpv << this->region << "\t"
                  << sv_it->first + 1 << "\t"
                  << ".\t"
                  << refbase << "\t"
                  << high_base << "\t"
                  << (short)snpq << "\t"
                  << filter << "\t"
                  //<< "P=" << pr_SNP_S_j_c_j << "\t"
                  << info_ss.str() << "\t"
                  << "GT:VR:RR:DP:GQ\t"
                  << genotype << ":"
                  << var_cov << ":"
                  << refbase_cov << ":"
                  << dp_cov << ":"
                  << ".\n";
        }
        else
            this->coverages->setP(sv_it->first, QUALSIDX_SNP_P, pr_SNP_S_j_c_j);
    }

    this->events->snps.erase(this->events->snps.begin(), sv_it);
}

void VcfWriter::printIndelBuffer(pos_t next_var_pos, const bed_coord_t &seg, std::ofstream &oindelv)
{
    char genotype[4];
    std::string filter, alt, ref;
    IndelMap::iterator iv;

    for (iv = this->events->indels.begin(); iv != this->events->indels.end() && iv->first < next_var_pos; ++iv)
    {
        coverage_t ar_cov = 0;
        auto ivg_it = iv->second.begin();
        IndelEvent &max_reads_indel = ivg_it->second;

        while (ivg_it != iv->second.end())
        {
            ar_cov += ivg_it->second.read_count;
            if (ivg_it->second.read_count > max_reads_indel.read_count)
                max_reads_indel = ivg_it->second;
            ++ivg_it;
        }

        pos_t var_start = max_reads_indel.var_start;
        if (!max_reads_indel.isdel)
            --var_start;
        coverage_t total_coverage = this->coverages->getCoverage(var_start, COVSIDX_INDEL_DP);
        double mpr_strand_dir = max_reads_indel.simpleStrandTest() ? 1.0 : 0.0;
        double mean_avnqs = max_reads_indel.getMeanAvnqs();
        double mpr_entropy = max_reads_indel.simpleLocalEntropy(*this->sequences);
        double PR_INDEL_j;

        if (this->opts->indeltest)
        {
            double mean_mapq = max_reads_indel.getMeanMapq();
            double mean_var_rate = max_reads_indel.getMeanVarRate();
            PR_INDEL_j = indel_logit_wgs(this->ill, mpr_entropy, mpr_strand_dir, mean_avnqs, mean_mapq, mean_var_rate, max_reads_indel.read_count, total_coverage);
        }
        else
            PR_INDEL_j = indel_logit(this->ill, mpr_entropy, mpr_strand_dir, mean_avnqs, max_reads_indel.read_count, total_coverage);

        qual_t qual = (1 - PR_INDEL_j < 0.0000000001) ? 100 : (double)round(-10.0 * log10(1.0 - PR_INDEL_j));
        //double p = (double)round(PR_INDEL_j * 10000.0) / 10000.0;
        double p = PR_INDEL_j;

        if (this->opts->capturebed != nullptr ||
            (max_reads_indel.var_start >= seg.first &&
             max_reads_indel.var_start <= seg.second))
        {
            if (max_reads_indel.isdel)
                this->coverages->addCoverageRange(var_start, COVSIDX_INDEL_VR, max_reads_indel.var_len, max_reads_indel.read_count);

            if (this->opts->allnonref ||
                qual == 0 ||
                p <= this->opts->INDEL_PR_CUTOFF ||
                (this->opts->SLX_SNP_HIGH_MAP_QUAL_MAX_COVERAGE != 0 &&
                 total_coverage > this->opts->SLX_SNP_HIGH_MAP_QUAL_MAX_COVERAGE))
            {
                continue;
            }

            filter.clear();
            double var_ratio = (double)max_reads_indel.read_count / total_coverage;

            if (max_reads_indel.read_count < this->opts->MIN_VAR_READS)
                addField(filter, "low_VariantReads");

            if (total_coverage < this->opts->MIN_DEPTH_COVERAGE)
                addField(filter, "low_coverage");
            else if (total_coverage > this->opts->SLX_SNP_HIGH_MAP_QUAL_MAX_COVERAGE)
                addField(filter, "high_coverage");

            if (var_ratio < this->opts->MIN_VAR_RATIO)
                addField(filter, "low_VariantRatio");
            if (this->opts->STRAND_DIR_FILTER && !max_reads_indel.simpleStrandTest())
                addField(filter, "single_strand");
            if (max_reads_indel.getNearReadEndRatio() > this->opts->MAX_NEAR_READ_END_RATIO)
                addField(filter, "read_end_ratio");

            *genotype = '\0';
            coverage_t n_tmp = 0;
            coverage_t total_ref_cov = (max_reads_indel.isdel)
                                           ? this->coverages->getCoverage(var_start, COVSIDX_INDEL_RR)
                                           : this->coverages->getCoverage(var_start, COVSIDX_INDEL_RR_INS) - ar_cov;

            if (max_reads_indel.read_count == 0)
            {
                if (total_coverage < this->opts->INDEL_DEPTH_CUTOFF)
                    strncpy(genotype, "./.", 4);
                else
                    n_tmp = total_coverage;
            }
            else
            {
                if (total_coverage < this->opts->INDEL_DEPTH_CUTOFF)
                    strncpy(genotype, "1/.", 4);
                else
                    n_tmp = max_reads_indel.read_count + total_ref_cov;
            }

            if (*genotype == '\0')
            {
                double ratio = (double)max_reads_indel.read_count / n_tmp;
                if (ratio >= this->opts->INDEL_HOM_VAR_CUTOFF)
                    strncpy(genotype, "1/1", 4);
                else if (ratio >= this->opts->INDEL_HET_CUTOFF)
                    strncpy(genotype, "1/0", 4);
                else
                {
                    ratio = (double)total_ref_cov / n_tmp;
                    if (ratio >= this->opts->INDEL_HOM_VAR_CUTOFF)
                        strncpy(genotype, "0/0", 4);
                    else if (ratio >= this->opts->INDEL_HET_CUTOFF)
                        strncpy(genotype, "0/.", 4);
                    else
                        strncpy(genotype, "./.", 4);
                }
            }

            ref.clear();
            alt.clear();

            if (max_reads_indel.isdel)
            {
                ref.append(this->sequences->getRefSeqRange(max_reads_indel.var_start - 1, max_reads_indel.var_len + 1));
                alt.push_back(this->sequences->seq[max_reads_indel.var_start - 1]);
            }
            else
            {
                char refbase = this->sequences->seq[max_reads_indel.var_start - 1];
                ref.push_back(refbase);
                alt.push_back(refbase);
                alt.append(max_reads_indel.seq);
            }

            if (p < 0.5)
                addField(filter, "low_qual");
            else if (filter.empty())
                filter = "PASS";

            if (this->opts->gvcftest && this->events->last_call_indel < max_reads_indel.var_start - 1)
            {
                outputGvcf(max_reads_indel.var_start - 1,
                           this->events->last_call_indel,
                           seg.second,
                           oindelv,
                           QUALSIDX_INDEL_P,
                           COVSIDX_INDEL_VR,
                           COVSIDX_INDEL_RR,
                           COVSIDX_INDEL_DP);
            }
            this->events->last_call_indel = max_reads_indel.var_start;

            oindelv << this->region << "\t"
                    << max_reads_indel.var_start << "\t"
                    << ".\t"
                    << ref << "\t"
                    << alt << "\t"
                    << (short)qual << "\t"
                    << filter << "\t"
                    << "P=" << p << "\t"
                    << "GT:VR:RR:DP:GQ\t"
                    << genotype << ":"
                    << max_reads_indel.read_count << ":"
                    << total_ref_cov << ":"
                    << total_coverage << ":"
                    << ".\n";
        }
        else
            this->coverages->setP(var_start, QUALSIDX_INDEL_P, p);
    }

    this->events->indels.erase(this->events->indels.begin(), iv);
}
