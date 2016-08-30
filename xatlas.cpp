#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include "htslib/faidx.h"
#include "htslib/sam.h"
/*#include "samtools/bam.h"*/
#include <map>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <bitset>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <signal.h>
#include <string>
#include <cstring>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <execinfo.h>
#include <sys/resource.h>
#include <utility>
#ifdef _OPENMP
#include <omp.h>
#endif
#define MEM_USAGE(X)                                                                                                   \
    {                                                                                                                  \
        char rar[MAX_READ_NAME];                                                                                       \
        snprintf(rar, MAX_READ_NAME, "echo 'MEM: %d ' ", X);                                                           \
        system(rar);                                                                                                   \
        snprintf(rar, MAX_READ_NAME,                                                                                   \
                 "ps -o cmd,pid,ppid,stat,bsdtime,stime,pcpu,size,rss,vsize,cputime,etime -p %d | tail -n2;",          \
                 getpid());                                                                                            \
        system(rar);                                                                                                   \
    }
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#ifndef __USE_GNU
#define __USE_GNU
#endif
#ifndef BITS
#define BITS "testing version"
#endif
#ifndef MAX_READ_NAME
#define MAX_READ_NAME 100
#endif
__asm__(".align\t8\n"
        ".LCCC2:\n\t"
        ".long\t0\n\t"
        ".long\t1072693248\n\t");
typedef std::vector<std::pair<uint32_t, uint32_t> > LIST;
typedef std::map<std::string, LIST> LISTOREGIONS;
void tokenise(std::vector<std::string>& t, std::string const& l, char s) {
    using namespace std;
    size_t pos = 0;
    size_t lst = 0;
    while (pos != std::string::npos) {
        pos = l.find(s, pos + 1);
        size_t g = lst == 0 ? 0 : 1;
        string a = l.substr(lst + g, pos - lst - g);
        t.push_back(a);
        lst = pos;
    }
}
static char const* VERSION = "v0.0.2-rc4";
static char const* MY_NAME = "xatlas";
static char const* OPTIONS = "Required arguments:\n"
                             "\n\t--ref          fasta format reference"
                             "\n\t--bam          sorted, indexed bam file"
                             "\n\t--samplename   sample name as appears in vcf"
                             "\n\t--prefix       output file prefix\n"
                             "\nOptions:\n"
                             "\n\t--capturebed   bed file of regions to inspect - only of use when"
                             "\n\t               combined with --gvcftest or --novar"
                             "\n\t--gvcf         gvcf [testing]"
                             "\n\t--novar        give metrics for uncalled sites"
                             "\n\t--allnonref    turn off filtering - i.e. call all non-ref sites"
                             "\n\t               (unlike no --novar will give allele but bit silly...)"
                             "\n\t--maxcov       maximum coverage events to report (8000)"
                             "\n\t--xindel       6-parameter model trained on HiSeqX NA12878 data"
                             "\n\t--blockabslim"
                             "\n\t--blockrellim"
                             "\n\t--nonreffrac"
                             "\n\t--nonreffraccutoff"
                             "\n\t--minindelmapq"
                             "\n\t--minsnpmapq\n";
namespace utils {
void touchfile(const std::string& file) {
    FILE* fp = fopen(file.data(), "ab+");
    fclose(fp);
}
bool isregfile(const char* fn) {
    struct stat test;
    if (stat(fn, &test) != 0) {
        return false;
    }
    return S_ISREG(test.st_mode);
}
bool isregfile(std::string a) { return isregfile(a.data()); }
void removefile(std::string const& tmp) {
    struct stat sdir;
    if (stat(tmp.data(), &sdir) != 0) {
        std::cout << "problem removing file " << tmp << "\n";
        if (remove(tmp.data()) != 0)
            std::cout << "problem removing file " << tmp << "\n";
    }
}
}
namespace {
typedef uint16_t COVERAGE_t;
typedef uint8_t QUAL_t;
typedef uint16_t BIG_QUAL_t;
typedef std::string KEY_t;
typedef std::string S_SEQ_t;
typedef char* CP_SEQ_t;
typedef uint32_t REF_ID_t;
typedef uint32_t POS_t;
typedef POS_t LENGTH_t;
typedef POS_t READ_POS_t;
typedef POS_t REF_POS_t;
typedef POS_t OFFSET_t;
typedef OFFSET_t READ_OFFSET_t;
typedef OFFSET_t READ_OFFSET_t;
typedef char BASE_t;
typedef uint32_t GENOTYPE;
size_t const GENOTYPE_LENGTH = 4;
typedef std::string RG_t;
typedef std::map<REF_POS_t, double> TMP_P_HOLDER;
static REF_POS_t last_call_indel;
static REF_POS_t last_call_snp;
}
namespace opts {
bool const compatibility = true;
bool snp_file = false;
bool checkpoint = true;
bool force = false;
bool depth = false;
bool rate = false;
bool indeltest = false;
bool allindels = false;
bool allnonref = false;
bool gvcftest = false;
bool gvcftestdebug = false;
bool run_in_fork = false;
bool silly_scavenge = true;
bool ignore_low_qual = true;
COVERAGE_t const MIN_VAR_READS = 2;
COVERAGE_t const MIN_TOTAL_DEPTH = 5;
bool const STRAND_DIR_FILTER = false;
QUAL_t const MIN_MAP_QUAL = 0;
BIG_QUAL_t const INDEL_QUAL_CUTOFF = 2;
COVERAGE_t const MIN_DEPTH_COVERAGE = 5;
double const MIN_VAR_RATIO = 0.06;
double const SNP_MIN_PR = 0.109;
double const INDEL_PR_CUTOFF = SNP_MIN_PR;
double const MAX_NEAR_READ_END_RATIO = 0.8;
COVERAGE_t const INDEL_DEPTH_CUTOFF = 5;
double const INDEL_HOM_VAR_CUTOFF = 0.6;
double const INDEL_HET_CUTOFF = 0.06;
double const SNP_HET_MIN = 0.1;
double const SNP_HET_MAX = 0.9;
double const SNP_STRAND_RAIO_CUTOFF = 0.01;
COVERAGE_t const SNP_STRAND_TEST_COV_CUTOFF = 16;
COVERAGE_t const SNP_MIN_COV = 2;
double const READ_LEVEL_Z_CUTOFF = -1.3;
unsigned block_abs_lim = 10;
double block_rel_lim = 0.3;
double block_rel_min = 1.0;
double max_alt_frac = 0.50;
COVERAGE_t VRCUTOFF = 0;
char block_label[25] = "min";
char const* known_indels = 0;
char const* capturebed = 0;
COVERAGE_t SLX_SNP_HIGH_MAP_QUAL_MAX_COVERAGE = 8000;
}
namespace {
typedef std::pair<double, double> mean_and_sd_vals;
typedef struct mean_and_sdx {
    unsigned long long k_;
    double mean_;
    double sd_;
    double min_;
    double max_;
    double block_rel_lim;
    unsigned block_abs_lim;
    unsigned block_rel_min;
} mean_and_sdx;
int breakblock2(mean_and_sdx* x, double v) {
    if (x->k_ == 0)
        return 0;
    int r = 0;
    __asm__ __volatile__("xorl\t%[ebx], %[ebx]\n\t"
                         ".cfi_def_cfa_offset 16\n\t"
                         ".cfi_offset 3, -16\n\t"
                         "movl\t52(%[rdi]), %%eax\n\t"
                         "movsd\t24(%[rdi]), %%xmm1\n\t"
                         "cvtsi2sdq\t%%rax, %%xmm2\n\t"
                         "ucomisd\t%%xmm2, %%xmm1\n\t"
                         "jb\t.LLLL1112%=\n\t"
                         "movsd\t40(%[rdi]), %%xmm2\n\t"
                         "movl\t$1, %[ebx]\n\t"
                         "mulsd\t%%xmm1, %%xmm2\n\t"
                         "ucomisd\t%%xmm2, %[xmm0]\n\t"
                         "ja\t.LLLL6%=\n\t"
                         ".LLLL1112%=:\n\t"
                         "movl\t48(%[rdi]), %%eax\n\t"
                         "xorl\t%[ebx], %[ebx]\n\t"
                         "cvtsi2sdq\t%%rax, %%xmm2\n\t"
                         "addsd\t%%xmm1, %%xmm2\n\t"
                         "ucomisd\t%%xmm2, %[xmm0]\n\t"
                         "seta\t%%bl\n\t"
                         ".LLLL6%=:\n\t"
                         : [ebx] "=b"(r)
                         : [xmm0] "x"(v), [rdi] "D"(x)
                         : "%rax", "%xmm1", "%xmm2");
    return r;
}
int breakblock3(mean_and_sdx* x, double v) {
    int r = 0;
    __asm__ __volatile__("xorl\t%k[eax], %k[eax]\n\t"
                         "cmpq\t$0, (%[rdi])\n\t"
                         "je\t.LLL45%=\n\t"
                         "movl\t52(%[rdi]), %k[eax]\n\t"
                         "movsd\t24(%[rdi]), %%xmm1\n\t"
                         "cvtsi2sdq\t%q[eax], %%xmm2\n\t"
                         "ucomisd\t%%xmm2, %%xmm1\n\t"
                         "jb\t.LLL46%=\n\t"
                         "movsd\t40(%[rdi]), %%xmm3\n\t"
                         "movl\t$1, %k[eax]\n\t"
                         "movsd\t.LCCC2(%%rip), %%xmm4\n\t"
                         "movapd\t%%xmm3, %%xmm2\n\t"
                         "addsd\t%%xmm4, %%xmm2\n\t"
                         "mulsd\t%%xmm1, %%xmm2\n\t"
                         "ucomisd\t%%xmm2, %[xmm0]\n\t"
                         "ja\t.LLL45%=\n\t"
                         "movapd\t%%xmm4, %%xmm2\n\t"
                         "subsd\t%%xmm3, %%xmm2\n\t"
                         "mulsd\t%%xmm1, %%xmm2\n\t"
                         "ucomisd\t%[xmm0], %%xmm2\n\t"
                         "ja\t.LLL45%=\n\t"
                         ".LLL46%=:\n\t"
                         "subsd\t%%xmm1, %[xmm0]\n\t"
                         "cvttsd2si\t%[xmm0], %k[eax]\n\t"
                         "cltd\n\t"
                         "xorl\t%%edx, %k[eax]\n\t"
                         "subl\t%%edx, %k[eax]\n\t"
                         "cmpl\t48(%[rdi]), %k[eax]\n\t"
                         "setnb\t%b[eax]\n\t"
                         "movzbl\t%b[eax], %k[eax]\n\t"
                         ".LLL45%=:\n\t"
                         : [eax] "=q"(r)
                         : [xmm0] "x"(v), [rdi] "D"(x)
                         : "%rdx", "%xmm1", "%xmm2", "%xmm3", "%xmm4");
    return r;
}
void add_val(mean_and_sdx* x, double v) {
    __asm__ __volatile__("movq\t(%[rdi]), %%rax\n\t"
                         "testq\t%%rax, %%rax\n\t"
                         "je\t.LL11125%=\n\t"
                         "movsd\t24(%[rdi]), %%xmm1\n\t"
                         "ucomisd\t%[xmm0], %%xmm1\n\t"
                         "ja\t.LL11126%=\n\t"
                         "ucomisd\t32(%[rdi]), %[xmm0]\n\t"
                         "jbe\t.LL1115%=\n\t"
                         "movsd\t%[xmm0], 32(%[rdi])\n\t"
                         ".LL1115%=:\n\t"
                         "addq\t$1, %%rax\n\t"
                         "movsd\t8(%[rdi]), %%xmm3\n\t"
                         "movapd\t%[xmm0], %%xmm2\n\t"
                         "testq\t%%rax, %%rax\n\t"
                         "movq\t%%rax, (%[rdi])\n\t"
                         "subsd\t%%xmm3, %%xmm2\n\t"
                         "js\t.LL1119%=\n\t"
                         "cvtsi2sdq\t%%rax, %%xmm1\n\t"
                         ".LL11120%=:\n\t"
                         "movapd\t%%xmm2, %%xmm4\n\t"
                         "divsd\t%%xmm1, %%xmm4\n\t"
                         "movapd\t%%xmm4, %%xmm1\n\t"
                         "addsd\t%%xmm3, %%xmm1\n\t"
                         "subsd\t%%xmm1, %[xmm0]\n\t"
                         "movsd\t%%xmm1, 8(%[rdi])\n\t"
                         "mulsd\t%%xmm2, %[xmm0]\n\t"
                         "addsd\t16(%[rdi]), %[xmm0]\n\t"
                         "movsd\t%[xmm0], 16(%[rdi])\n\t"
                         "jmp\t.LL111444%=\n\t"
                         ".p2align 4,,10\n\t"
                         ".p2align 3\n\t"
                         ".LL11126%=:\n\t"
                         "movsd\t%[xmm0], 24(%[rdi])\n\t"
                         "jmp\t.LL1115%=\n\t"
                         ".p2align 4,,10\n\t"
                         ".p2align 3\n\t"
                         ".LL11125%=:\n\t"
                         "movsd\t%[xmm0], 24(%[rdi])\n\t"
                         "movsd\t%[xmm0], 32(%[rdi])\n\t"
                         "jmp\t.LL1115%=\n\t"
                         ".p2align 4,,10\n\t"
                         ".p2align 3\n\t"
                         ".LL1119%=:\n\t"
                         "movq\t%%rax, %%rdx\n\t"
                         "andl\t$1, %%eax\n\t"
                         "shrq\t%%rdx\n\t"
                         "orq\t%%rax, %%rdx\n\t"
                         "cvtsi2sdq\t%%rdx, %%xmm1\n\t"
                         "addsd\t%%xmm1, %%xmm1\n\t"
                         "jmp\t.LL11120%=\n\t"
                         ".LL111444%=:\n\t"
                         :
                         : [xmm0] "x"(v), [rdi] "S"(x)
                         : "%rax", "%rdx", "%xmm1", "%xmm2", "%xmm3", "%xmm4");
    return;
}
#define getmin_(X) ((X)->min_)
#define getmax_(X) ((X)->max_)
#define getmean_(X) ((X)->mean_)
#define getk_(X) ((X)->k_)
inline void reset_blocker(mean_and_sdx* x) {
    memset(x, 0, sizeof(mean_and_sdx));
    x->block_rel_lim = opts::block_rel_lim;
    x->block_abs_lim = opts::block_abs_lim;
    x->block_rel_min = opts::block_rel_min;
}
inline double getsd_(mean_and_sdx* x) {
    double t = x->sd_;
    if (x->k_ >= 2)
        t = pow(t / (x->k_ - 1), 0.5);
    return t;
}
typedef std::map<std::string, double> KNOWN_INDELS;
KNOWN_INDELS get_known(std::string& region) {
    using namespace std;
    std::map<string, double> list;
    cout << "grab known\n";
    ifstream pvcf(opts::known_indels);
    if (!pvcf)
        cout << "couldn't open dbsnp file?!?", exit(1);
    string line;
    unsigned z = 0;
    while (getline(pvcf, line)) {
        if (line[0] == '#')
            continue;
        vector<string> t;
        tokenise(t, line, '\t');
        if (t[0] != region)
            continue;
        size_t f = t[4].find(",");
        if (f != string::npos)
            continue;
        string type = t[1];
        if (t[3] < t[4]) {
            type += "I";
            type += t[4].substr(1, (int)(t[4].length()));
        } else {
            type += "D";
            char len[5];
            snprintf(len, 5, "%d", (int)(t[3].length() - 1));
            type += len;
        }
        f = t[7].find("GMAF=");
        if (f != string::npos) {
            size_t n = t[7].find(";", f + 1);
            list.insert(make_pair(type, atof(t[7].substr(f + 5, n - f - 5).data())));
        } else {
            list.insert(make_pair(type, 0.0001));
        }
        ++z;
    }
    cout << "grabbed " << z << " from " << region << "\n";
    return list;
}
}
namespace cuidado {
static uint32_t const MAGIC = 255 * 1024 * 1024;
size_t const MAX_REF_NAME_LENGTH = 10;
inline void switch_it(char x, std::string& y, size_t z) {
    switch (x) {
    case 'a':
        y[z] = 't';
        break;
    case 'c':
        y[z] = 'g';
        break;
    case 'g':
        y[z] = 'c';
        break;
    case 't':
        y[z] = 'a';
        break;
    case 'n':
        y[z] = 'n';
        break;
    case 'A':
        y[z] = 'T';
        break;
    case 'C':
        y[z] = 'G';
        break;
    case 'G':
        y[z] = 'C';
        break;
    case 'T':
        y[z] = 'A';
        break;
    case 'N':
        y[z] = 'N';
        break;
    default:
        assert(0);
        break;
    }
}
}
namespace ill {
double const SLX_SNP_intercept = -9.088;
double const SLX_SNP_quality_score = 0.162;
double const SLX_SNP_NQS_pass = 1.645;
double const SLX_SNP_swap = 0.0;
double const SLX_SNP_rel_pos = 2.349;
double const SLX_SNP_MAX_SUB = 5.0;
double const SLX_SNP_MAX_INDEL = 5.0;
COVERAGE_t const SLX_INDEL_LOW_MAP_QUAL_MAX_COVERAGE = ((COVERAGE_t)-1) >> 2;
COVERAGE_t const SLX_SNP_HIGH_MAP_QUAL_MAX_COVERAGE = 8000;
COVERAGE_t const SLX_SNP_MIN_COVERAGE = 6;
double const prior_err_c = 0.1;
double const prior_snp_c = 0.9;
double const SLX_SNP_CUTOFF = 0.95;
double const SLX_err_prior_arr[] = {0.991, 0.004, 0.002, 0.001, 0.001, 0.001, 0.000002, 0.000005, 0.0000001, 0.000001};
double const SLX_snp_prior_arr[] = {0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001,
                                    0.0000001, 0.0000001, 0.0000001, 0.027,     0.973};
inline double snp_logit(double q, double s, double n, double r) {
    double tmp = SLX_SNP_intercept + (q * SLX_SNP_quality_score) + (s * SLX_SNP_NQS_pass) + (n * SLX_SNP_swap) +
                 (r * SLX_SNP_rel_pos);
    return (exp(tmp) / (1 + exp(tmp)));
}
double const SLX_INDEL_intercept = -20.5000;
double const SLX_INDEL_local_entropy = 3.39065;
double const SLX_INDEL_strand_dir = 3.02573;
double const SLX_INDEL_norm_var_sq = 0.32695;
double const SLX_INDEL_mean_avenqs = 0.37184;
double const WGS_INDEL_intercept = -14.95535266;
double const WGS_INDEL_simple_local_entropy = 0.16228492;
double const WGS_INDEL_strand_dir = 4.84302146;
double const WGS_INDEL_norm_var_square = 0.02302639;
double const WGS_INDEL_mean_ave_nqs = 0.18140485;
double const WGS_INDEL_mean_map_qual = 0.10577339;
double const WGS_INDEL_mean_var_rate = -62.03643929;
inline double indel_logit_wgs(double local_ref_simple_entropy, double strand_dir, double mean_avnqs, double mean_mapq,
                              double mean_var_rate, COVERAGE_t read_count, COVERAGE_t total_depth) {
    double tmp = WGS_INDEL_intercept + (local_ref_simple_entropy * WGS_INDEL_simple_local_entropy) +
                 (strand_dir * WGS_INDEL_strand_dir) + (mean_avnqs * WGS_INDEL_mean_ave_nqs) +
                 (WGS_INDEL_norm_var_square * (((double)pow(read_count, 2.0)) / total_depth)) +
                 (WGS_INDEL_mean_map_qual * mean_mapq) + (WGS_INDEL_mean_var_rate * mean_var_rate);
    if (tmp > 700)
        return 1.0;
    return (exp(tmp) / (1 + exp(tmp)));
}
inline double indel_logit(double local_ref_simple_entropy, double strand_dir, double mean_avnqs, COVERAGE_t read_count,
                          COVERAGE_t total_depth) {
    double tmp = SLX_INDEL_intercept + (local_ref_simple_entropy * SLX_INDEL_local_entropy) +
                 (strand_dir * SLX_INDEL_strand_dir) + (mean_avnqs * SLX_INDEL_mean_avenqs) +
                 ((((double)pow(read_count, 2.0)) / total_depth) * SLX_INDEL_norm_var_sq);
    if (tmp > 700)
        return 1.0;
    return (exp(tmp) / (1 + exp(tmp)));
}
}
namespace VCF {
struct VAR {
    S_SEQ_t scf_;
    std::string tmpinfo_;
    REF_POS_t pos_;
    S_SEQ_t ref_;
    S_SEQ_t alt_;
    S_SEQ_t filters_;
    QUAL_t qual_;
    GENOTYPE GT_[GENOTYPE_LENGTH];
    uint32_t VR_;
    uint32_t AR_;
    uint32_t RR_;
    uint32_t DP_;
    double pr_;
    bool used;
    VAR(S_SEQ_t s, REF_POS_t p, S_SEQ_t r, S_SEQ_t a, S_SEQ_t f, QUAL_t q, char* g, uint16_t vr, uint16_t ar,
        uint16_t rr, uint16_t dp, double pr, bool u, std::string i)
        : scf_(s), tmpinfo_(i), pos_(p), ref_(r), alt_(a), filters_(f), qual_(q), VR_(vr), AR_(ar), RR_(rr), DP_(dp),
          pr_(pr), used(u) {
        memset(GT_, 0, sizeof(GENOTYPE) * GENOTYPE_LENGTH);
        memcpy(GT_, g, GENOTYPE_LENGTH - 1);
    }
};
std::ostream& operator<<(std::ostream& os, const VAR& v) {
    os << v.scf_ << "\t" << v.pos_ << "\t"
       << (v.ref_.length() > v.alt_.length() ? "DEL" : v.ref_.length() < v.alt_.length() ? "INS" : "SNP") << "\t"
       << v.ref_ << "\t" << v.alt_ << "\tQ=" << (short)v.qual_ << "\t" << v.filters_ << "\tP=" << v.pr_
       << "\tVR=" << v.VR_ << "\tAR=" << v.AR_ << "\tRR=" << v.RR_ << "\tDP=" << v.DP_ << "\t" << (char*)v.GT_
       << "\tC=" << v.used << "\t" << v.tmpinfo_ << "\n";
    return os;
}
}
namespace {
using namespace std;
struct SNP_INFO {
    std::vector<uint64_t> readids_;
    std::map<RG_t, COVERAGE_t> readgroups_;
    SNP_INFO() : readids_(), readgroups_() {}
};
typedef std::map<BASE_t, SNP_INFO> SNP_INFOS;
class STUFF {
    bam_hdr_t* header;
    hts_idx_t* idx;
    STUFF();

  public:
    samFile* in;
    bam1_t* alignment;
    hts_itr_t* iter;
    STUFF(char const* bf) {
        in = sam_open_format(bf, "rb", 0);
        assert(in != 0);
        header = sam_hdr_read(in);
        assert(header != 0);
        idx = sam_index_load(in, bf);
        assert(idx != 0);
        alignment = bam_init1();
    }
    bam_hdr_t* get_header() const { return header; }
    hts_itr_t* get_iter(char const* what) {
        iter = sam_itr_querys(idx, header, what);
        return iter;
    }
    ~STUFF() {
        hts_itr_destroy(iter);
        hts_idx_destroy(idx);
        bam_hdr_destroy(header);
        bam_destroy1(alignment);
        sam_close(in);
    }
};
class SCF_SEQ {
    faidx_t* fai_;
    char* scf_seq_;
    char const* scf_name_;
    char const* scf_file_;
    uint32_t scf_len_;
    char get_nt(REF_POS_t i) const { return scf_seq_[i - 1]; }

  public:
    SCF_SEQ(char const* refseq, char const* region) : scf_name_(region), scf_file_(refseq) {
        int l;
        fai_ = fai_load(refseq);
        if (fai_ == 0)
            exit(1);
        scf_seq_ = fai_fetch(fai_, region, &l);
        scf_len_ = l;
        char* tmp = new char[l + 1];
        tmp[l] = '\0';
        memcpy(tmp, scf_seq_, l * sizeof(char));
        free(scf_seq_);
        scf_seq_ = tmp;
    }
    ~SCF_SEQ() {
        delete[] scf_seq_;
        fai_destroy(fai_);
    }
    BASE_t get_ref_base(READ_OFFSET_t n) const { return (get_nt(n)); }
    S_SEQ_t get_ref_seq_range(READ_OFFSET_t n, LENGTH_t l) const {
        std::string tmp;
        for (size_t i = n; i < n + l; ++i)
            tmp += (BASE_t)::toupper(get_nt(i));
        return tmp;
    }
    uint32_t get_length() const { return scf_len_; }
};
class COV_COUNTER {
    union COVS {
        uint64_t all_;
        uint16_t depths_[4];
    };
    union FCALLS {
        uint64_t all_;
        uint16_t depths_[4];
    };
    union QUALS {
        float quals_[2];
    };
    COVS* const covs_;
    FCALLS* const fcalls_;
    QUALS* const quals_;
    std::bitset<cuidado::MAGIC>& seen_;
    uint32_t len_;

  public:
    enum INDEX { AID_DEP = 0, AID_REF_DEPTH = 1, AS_COV = 2, AS_DEL = 3, DUMMY };
    enum FORCE_INDEX { AS_VR = 0, AS_AR = 1, AID_VR = 2, AID_AR = 3 };
    enum QUALS_INDEX { SNVP = 0, INDELP = 1 };
    COV_COUNTER(uint32_t l)
        : covs_(new COVS[l + 1]), fcalls_(new FCALLS[l + 1]), quals_(new QUALS[l + 1]),
          seen_(*(new std::bitset<cuidado::MAGIC>)), len_(l) {
        memset(quals_, 0, (l + 1) * sizeof(QUALS));
        memset(covs_, 0, (l + 1) * sizeof(COVS));
        memset(fcalls_, 0, (l + 1) * sizeof(FCALLS));
    }
    void set_snp_p(REF_POS_t p, float q) const { quals_[p].quals_[SNVP] = q; }
    void set_indel_p(REF_POS_t p, float q) const { quals_[p].quals_[INDELP] = q; }
    float get_p(REF_POS_t p, QUALS_INDEX g) const { return quals_[p].quals_[g]; }
    float get_snp_p(REF_POS_t p) const { return quals_[p].quals_[SNVP]; }
    float get_indel_p(REF_POS_t p) const { return quals_[p].quals_[INDELP]; }
    ~COV_COUNTER() {
        delete[] covs_;
        delete[] fcalls_;
        delete quals_, delete &seen_;
    }
    void add_coverage(INDEX i, uint32_t p) {
        if (covs_[p].depths_[i] < ill::SLX_INDEL_LOW_MAP_QUAL_MAX_COVERAGE)
            covs_[p].depths_[i]++;
    }
    void subtract_coverage(INDEX i, uint32_t p) { covs_[p].depths_[i]--; }
    void set_coverage(INDEX i, uint32_t p, uint16_t c) {
        covs_[p].depths_[i] =
            (c >= ill::SLX_INDEL_LOW_MAP_QUAL_MAX_COVERAGE ? ill::SLX_INDEL_LOW_MAP_QUAL_MAX_COVERAGE : c);
    }
    uint16_t get_coverage(INDEX i, uint32_t p) const { return covs_[p].depths_[i]; }
    bool has_depth(uint32_t p) const { return covs_[p].all_; }
    uint16_t is_toxic(INDEX i, uint32_t p) const {
        return covs_[p].depths_[i] >= ill::SLX_INDEL_LOW_MAP_QUAL_MAX_COVERAGE;
    }
    bool seen(uint32_t i) {
        bool tmp = seen_[i];
        if (!tmp)
            seen_[i] = true;
        return tmp;
    }
    void add_force_call(FORCE_INDEX i, uint32_t p) const {
        if (fcalls_[p].depths_[i] < (((uint16_t)-1) >> 4))
            fcalls_[p].depths_[i]++;
    }
    void set_force_call(FORCE_INDEX i, REF_POS_t p, COVERAGE_t c) const {
        fcalls_[p].depths_[i] = c > (((uint16_t)-1) >> 4) ? ((uint16_t)-1) : (uint16_t)c;
    }
    uint16_t get_force_call(FORCE_INDEX i, uint32_t p) const { return fcalls_[p].depths_[i]; }
    bool has_force_call(uint32_t p) const { return fcalls_[p].all_; }
};
bool indel_sorter(std::pair<S_SEQ_t, std::string> const& a, std::pair<S_SEQ_t, std::string> const& b) {
    return (a.second < b.second);
}
uint32_t longest_homopolymer_run(char const* refenv) {
    char last_base = '\0';
    uint32_t len = strlen(refenv) - 1;
    uint32_t* counts = new uint32_t[len];
    memset(counts, 0, (strlen(refenv) - 1) * sizeof(uint32_t));
    bool previous_char_same = false;
    uint32_t repeat_pos = 0;
    for (uint32_t i = 0; i < strlen(refenv); ++i) {
        char current_bast = ::toupper(refenv[i]);
        if (current_bast == last_base) {
            if (!previous_char_same) {
                repeat_pos = i - 1;
                counts[repeat_pos] = 2;
            } else {
                counts[repeat_pos]++;
            }
            previous_char_same = true;
        } else
            previous_char_same = false;
        last_base = current_bast;
    }
    uint32_t max = 1;
    for (uint32_t i = 0; i < len; ++i)
        if (counts[i] > max)
            max = counts[i];
    delete[] counts;
    return max;
}
inline int32_t cigar_op(uint32_t*& cit) { return (*cit & BAM_CIGAR_MASK); }
inline int32_t cigar_len(uint32_t*& cit) { return ((*cit & ~(BAM_CIGAR_MASK)) >> BAM_CIGAR_SHIFT); }
inline void cigar2asci(bam1_t*& b) {
    uint32_t* cit = bam_get_cigar(b);
    std::cout << "\ncigar2asci CIGAR= ";
    for (size_t i = 0; i < b->core.n_cigar; ++i, ++cit) {
        int32_t c_op = cigar_op(cit);
        REF_POS_t c_len = cigar_len(cit);
        switch (c_op) {
        case BAM_CMATCH:
            std::cout << c_len << "M";
            break;
        case BAM_CINS:
            std::cout << c_len << "I";
            break;
        case BAM_CDEL:
            std::cout << c_len << "D";
            break;
        case BAM_CSOFT_CLIP:
            std::cout << c_len << "S";
            break;
        default:
            assert(0);
            break;
        }
    }
    std::cout << "\n";
}
inline BASE_t get_read_base(bam1_t*& b, READ_OFFSET_t n) { return seq_nt16_str[bam_seqi(bam_get_seq(b), n)]; }
inline REF_POS_t get_read_ref_pos(bam1_t*& b) { return b->core.pos; }
inline REF_ID_t get_read_ref_id(bam1_t*& b) { return b->core.tid; }
inline S_SEQ_t get_read_seq_range(bam1_t*& b, READ_OFFSET_t n, LENGTH_t l) {
    std::string tmp;
    for (size_t i = n; i < n + l; ++i)
        tmp += get_read_base(b, i);
    return tmp;
}
int32_t const NEAR_END = 3;
int const SREADENV = 6;
int const READENV = (SREADENV * 2) + 2;
inline bool as_nqs(uint32_t qpos, uint32_t dist, uint8_t* qquals, uint32_t qlen) {
    bool as_pstat = false;
    if (qlen > 10) {
        if (dist > qlen - 6) {
            as_pstat = true;
        } else if (qpos >= 5 && qlen - qpos >= 6 && qquals[qpos] >= 20) {
            bool all = true;
            for (int i = -5; i < 6; ++i) {
                if (i == 0)
                    continue;
                if (qquals[qpos + i] < 15)
                    all = false;
            }
            as_pstat = all;
        }
    }
    return as_pstat;
}
struct SNP_EVENT {
  private:
    SNP_EVENT();

  public:
    SNP_EVENT(BASE_t ab, BASE_t rb, QUAL_t q, READ_POS_t qp, READ_OFFSET_t d3, bool nq, char const* re, char const* rn,
              bool rs, LENGTH_t rl, int score, double srate, double grate, uint64_t ri, RG_t rg)
        : readid_(ri), rg_(rg), read_len_(rl), read_sub_rate_(srate), read_gap_rate_(grate), read_score_(score),
          qplace_(qp), dist3_(d3), qual_(q), allele_base_(ab), ref_base_(rb), nqs_(nq), read_strand_(rs), type_(blah) {
        memset(read_readenv_, 0, READENV * sizeof(BASE_t));
        strcpy(read_readenv_, re);
        memset(read_name_, 0, MAX_READ_NAME * sizeof(BASE_t));
        strncpy(read_name_, rn, MAX_READ_NAME - 1);
    }
    char read_name_[MAX_READ_NAME];
    BASE_t read_readenv_[READENV];
    uint64_t readid_;
    RG_t rg_;
    LENGTH_t read_len_;
    double read_sub_rate_;
    double read_gap_rate_;
    int read_score_;
    READ_POS_t qplace_;
    READ_OFFSET_t dist3_;
    QUAL_t qual_;
    BASE_t allele_base_;
    BASE_t ref_base_;
    bool nqs_;
    bool read_strand_;
    enum { SNP_ = 0, SWAP_ = 1, MNP_ = 2, blah } type_;
    enum { BASE_A = 0, BASE_C = 1, BASE_G = 2, BASE_T = 3, BASE_N };
};
struct VAR_CORE {
    REF_ID_t ref_tid_;
    REF_POS_t var_start_;
    LENGTH_t var_len_;
    VAR_CORE() : ref_tid_(0), var_start_(0), var_len_(0) {}
    VAR_CORE(REF_ID_t r, REF_POS_t s, LENGTH_t l) : ref_tid_(r), var_start_(s), var_len_(l) {}
};
struct NEARBY_VARS : public VAR_CORE {
    enum TYPE_ { INS = 0, DEL, SNP, blah } type_;
    S_SEQ_t seq_;
    NEARBY_VARS() : VAR_CORE(0, 0, 0), type_(blah), seq_() {}
    NEARBY_VARS(REF_ID_t r, REF_POS_t s, TYPE_ t, BASE_t sq) : VAR_CORE(r, s, -1), type_(t), seq_() { seq_ += sq; }
    NEARBY_VARS(REF_ID_t r, REF_POS_t s, LENGTH_t l, TYPE_ t, S_SEQ_t sq) : VAR_CORE(r, s, l), type_(t), seq_(sq) {}
    NEARBY_VARS(REF_ID_t r, REF_POS_t s, LENGTH_t l, TYPE_ t) : VAR_CORE(r, s, l), type_(t), seq_("") {}
    KEY_t silly_key() const {
        std::stringstream os("");
        os << ref_tid_ << ":" << var_start_;
        if (type_ == NEARBY_VARS::INS)
            os << "I";
        else if (type_ == NEARBY_VARS::DEL)
            os << "D" << var_len_;
        os << seq_;
        return os.str();
    }
};
typedef std::vector<NEARBY_VARS> NEARBY_VAR_SITES;
typedef std::map<const REF_POS_t, SNP_EVENT> A_READS_SNP_LIST;
struct A_READS_INDEL : public VAR_CORE {
  private:
    A_READS_INDEL();

  public:
    A_READS_INDEL(REF_ID_t r, REF_POS_t p, std::string id, uint32_t ls, uint32_t ts, uint32_t ao, LENGTH_t l)
        : VAR_CORE(r, p, l), seq_(id), lsclip(ls), tsclip(ts), offset(ao) {}
    S_SEQ_t seq_;
    uint32_t lsclip;
    uint32_t tsclip;
    int32_t offset;
    KEY_t silly_key() const {
        std::stringstream os("");
        os << var_start_;
        if (seq_.empty())
            os << "D" << var_len_;
        else
            os << "I" << seq_;
        return os.str();
    }
};
std::ostream& operator<<(std::ostream& os, const SNP_EVENT& snp) {
    os << "[\"" << snp.allele_base_ << "\", " << (short)snp.qual_ << ", " << snp.qplace_ << ", " << snp.dist3_ << ", "
       << snp.nqs_ << ", \"" << snp.read_readenv_ << "\", "
       << (snp.type_ == SNP_EVENT::SNP_ ? "0" : snp.type_ == SNP_EVENT::SWAP_ ? "1" : snp.type_ == SNP_EVENT::MNP_
                                                                                          ? "0.5"
                                                                                          : "WTF") << "]";
    return os;
}
std::ostream& operator<<(std::ostream& os, const NEARBY_VARS& vs) {
    os << vs.var_start_;
    if (vs.type_ == NEARBY_VARS::INS)
        os << "I";
    else if (vs.type_ == NEARBY_VARS::DEL)
        os << "D" << vs.var_len_;
    os << vs.seq_;
    return os;
}
std::ostream& operator<<(std::ostream& os, const A_READS_INDEL& indel) {
    os << indel.silly_key();
    return os;
}
typedef std::vector<A_READS_INDEL> A_READS_INDEL_LIST;
typedef std::map<const REF_POS_t, COVERAGE_t> REF_COV_t;
typedef REF_COV_t DEPTH;
struct INDEL_EVENT : public VAR_CORE {
    std::vector<uint64_t> readids_;
    std::map<RG_t, COVERAGE_t> readgroups_;
    std::string tmpinfo_;
    KEY_t indel_;
    double var_rate_gap_and_mismatch_;
    uint16_t read_count_;
    uint16_t near_read_end_count_;
    BIG_QUAL_t ave_nbq_;
    BIG_QUAL_t map_qual_;
    QUAL_t max_map_qual_;
    bool has_calling_read_in_pos_strand_;
    bool has_calling_read_in_neg_strand_;
    bool required_include_;
    bool isdel_;
    INDEL_EVENT(REF_ID_t ref, REF_POS_t start, REF_POS_t len, bool pos, bool neg, uint16_t nearend, BIG_QUAL_t avenbq,
                QUAL_t mqual, KEY_t id, double vrate, bool il, uint64_t ri, RG_t rg)
        : VAR_CORE(ref, start, len), indel_(id), var_rate_gap_and_mismatch_(vrate), read_count_(1),
          near_read_end_count_(nearend), ave_nbq_(avenbq), map_qual_(mqual), max_map_qual_(mqual),
          has_calling_read_in_pos_strand_(pos), has_calling_read_in_neg_strand_(neg), required_include_(false),
          isdel_(il) {
        readids_.push_back(ri);
        readgroups_[rg] = 1;
    }
    INDEL_EVENT()
        : VAR_CORE(0, 0, 0), indel_(""), var_rate_gap_and_mismatch_(0), read_count_(0), near_read_end_count_(0),
          ave_nbq_(0), map_qual_(0), max_map_qual_(0), has_calling_read_in_pos_strand_(false),
          has_calling_read_in_neg_strand_(false), required_include_(false), isdel_(0) {}
    void add_read_indel_event(INDEL_EVENT const& iv) {
        ++read_count_;
        if (iv.has_calling_read_in_pos_strand_)
            has_calling_read_in_pos_strand_ = true;
        if (iv.has_calling_read_in_neg_strand_)
            has_calling_read_in_neg_strand_ = true;
        near_read_end_count_ += iv.near_read_end_count_;
        if (iv.map_qual_ > max_map_qual_)
            max_map_qual_ = iv.map_qual_;
        map_qual_ += iv.map_qual_;
        ave_nbq_ += iv.ave_nbq_;
        var_rate_gap_and_mismatch_ += iv.var_rate_gap_and_mismatch_;
        readids_.push_back(iv.readids_[0]);
        readgroups_[iv.readgroups_.begin()->first]++;
    }
    bool simple_strand_test() { return (has_calling_read_in_pos_strand_ && has_calling_read_in_neg_strand_); }
    double near_read_end_ratio() { return ((double)near_read_end_count_ / read_count_); }
    COVERAGE_t total_depth(DEPTH const& depth) { return depth.find(var_start_)->second; }
    double var_ratio(DEPTH const& depth) { return ((double)read_count_ / total_depth(depth)); }
    double mean_avnqs() { return (double)round(((double)ave_nbq_ / read_count_) * 100.0) / 100.0; }
    double mean_mapq() { return (double)round(((double)map_qual_ / read_count_) * 100.0) / 100.0; }
    double mean_var_rate() { return (double)round(((double)var_rate_gap_and_mismatch_ / read_count_) * 100.0) / 100.0; }
    double simple_local_entropy(SCF_SEQ const& sequences) {
        REF_POS_t lower = var_start_ < 11 ? 1 : isdel_ ? var_start_ - 9 : var_start_ - 10;
        REF_POS_t upper = isdel_ ? var_start_ + var_len_ : var_start_ + 1;
        upper += 11;
        S_SEQ_t seq = sequences.get_ref_seq_range(lower, upper - lower);
        if (seq.length() < var_len_)
            return 0.0;
        uint16_t* tmp = new uint16_t[2 * seq.length()];
        memset(tmp, 0, sizeof(uint16_t) * 2 * seq.length());
        char* p = new char[seq.length() + 1];
        memset(p, 0, (seq.length() + 1) * sizeof(char));
        char* q = new char[seq.length() + 1];
        memset(q, 0, (seq.length() + 1) * sizeof(char));
        REF_POS_t n = 0, i = 0;
        while (i + var_len_ < seq.length()) {
            bool found = false;
            strcpy(p, seq.data() + i);
            size_t j = 0;
            while (j < n) {
                bool match = true;
                strcpy(q, seq.data() + tmp[2 * j]);
                for (size_t k = 0; k < var_len_; k++) {
                    if (p[k] != q[k]) {
                        match = false;
                        break;
                    }
                }
                if (match) {
                    found = true;
                    tmp[j * 2 + 1] += 1;
                    break;
                }
                ++j;
            }
            if (!found) {
                tmp[n * 2] = i;
                tmp[n * 2 + 1] = 1;
                ++n;
            }
            ++i;
        }
        delete[] p;
        delete[] q;
        double e = 0.0;
        double s = 1.0 / (seq.length() - var_len_ + 1);
        for (size_t i2 = 0; i2 < n; ++i2) {
            double f = s * (double)tmp[i2 * 2 + 1];
            e -= f * log(f);
        }
        delete[] tmp;
        return (double)round((e)*1000.0) / 1000.0;
    }
};
struct NASTY_GLOBAL_SNP_HOLDER {
    typedef std::map<const REF_POS_t, std::vector<SNP_EVENT> > SNPS;
    typedef std::map<S_SEQ_t, INDEL_EVENT> INDEL_GROUP;
    typedef std::map<const REF_POS_t, INDEL_GROUP> INDELS;
    SNPS snps_;
    INDELS indel_;
    void add_snp(REF_POS_t p, SNP_EVENT const& s) { snps_[p].push_back(s); }
    void add_indel(REF_POS_t p, INDEL_EVENT const& i) {
        if (indel_.count(p) == 0 || indel_.find(p)->second.count(i.indel_) == 0)
            indel_[p].insert(make_pair(i.indel_, i));
        else
            indel_.find(p)->second.find(i.indel_)->second.add_read_indel_event(i);
    }
};
inline void output_novar(REF_POS_t rposit_first, REF_POS_t last_call_novar, LIST const& segs, uint32_t piece_i,
                         COV_COUNTER const& coverages, SCF_SEQ const& sequences, std::string const& region,
                         std::ofstream& osnpv, COV_COUNTER::QUALS_INDEX i_quals_, COV_COUNTER::FORCE_INDEX i_vr_,
                         COV_COUNTER::INDEX i_rr_, COV_COUNTER::INDEX i_dp_ = COV_COUNTER::DUMMY) {
    if (segs[piece_i].first > last_call_novar)
        last_call_novar = segs[piece_i].first;
    if (rposit_first != last_call_novar - 1) {
        for (REF_POS_t i = last_call_novar; i < rposit_first; ++i) {
            stringstream tmppb2;
            string silly = coverages.get_coverage(i_dp_, i) == 0 ? "no_Reads" : "no_VariantReads";
            BASE_t refbase = sequences.get_ref_base(i);
            tmppb2 << region << "\t" << i << "\t.\t" << refbase << "\t.\t" << coverages.get_p(i, i_quals_) << "\t"
                   << silly << "\t.\tGT:VR:RR:DP:GQ\t0/0:" << (short)coverages.get_force_call(i_vr_, i) << ":"
                   << coverages.get_coverage(i_rr_, i) << ":"
                   << (i_dp_ == COV_COUNTER::DUMMY
                           ? coverages.get_force_call(i_vr_, i) + coverages.get_coverage(i_rr_, i)
                           : (REF_POS_t)coverages.get_coverage(i_dp_, i)) << ":.\n";
            osnpv << tmppb2.str();
        }
    }
}
inline void output_gvcf(REF_POS_t rposit_first, REF_POS_t& last_call_gvcf, LIST const& segs, uint32_t piece_i,
                        COV_COUNTER const& coverages, SCF_SEQ const& sequences, std::string const& region,
                        std::string& output, COV_COUNTER::QUALS_INDEX i_quals_, COV_COUNTER::FORCE_INDEX i_vr_,
                        COV_COUNTER::INDEX i_rr_, COV_COUNTER::INDEX i_dp_ = COV_COUNTER::DUMMY) {
    bool lrate = rposit_first == sequences.get_length();
    if (rposit_first != last_call_gvcf - 1) {
        mean_and_sdx mas_q_;
        mean_and_sdx mas_vr_;
        mean_and_sdx mas_rr_;
        mean_and_sdx mas_dp_;
        reset_blocker(&mas_q_);
        reset_blocker(&mas_vr_);
        reset_blocker(&mas_rr_);
        reset_blocker(&mas_dp_);
        timespec nsstart, nsend;
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &nsstart);
        stringstream tmppb;
        for (REF_POS_t i = last_call_gvcf + 1, last_block = last_call_gvcf + 1; i < rposit_first; ++i) {
            if (opts::rate && lrate) {
                if (i % 10000 == 0) {
                    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &nsend);
                    long nsd = nsend.tv_nsec - nsstart.tv_nsec;
                    cerr << "final gvcf @ " << i << " (" << (double)10000000000000 / abs(nsd) << " pos/s)\r";
                    nsstart = nsend;
                }
            }
            REF_POS_t dp_cov = i_dp_ == COV_COUNTER::DUMMY
                                   ? coverages.get_force_call(i_vr_, i) + coverages.get_coverage(i_rr_, i)
                                   : (REF_POS_t)coverages.get_coverage(i_dp_, i);
            bool breakit = false;
            const char* reason = 0;
            string extra = "";
            if (i == segs[piece_i].second) {
                reason = "endofregion", breakit = true;
                last_call_gvcf = -1;
            }
            else if (i == rposit_first - 1)
                reason = "blocked", breakit = true;
            else if (breakblock3(&mas_dp_, (double)dp_cov))
                reason = "dp", breakit = true;
            else if (breakblock3(&mas_rr_, (double)coverages.get_coverage(i_rr_, i)))
                reason = "rr", breakit = true;
            else if (breakblock3(&mas_vr_, (double)coverages.get_force_call(i_vr_, i)))
                reason = "vr", breakit = true;
            else if (breakblock3(&mas_q_, (double)coverages.get_p(i, i_quals_)))
                reason = "QUAL", breakit = true;
            else if (getmin_(&mas_dp_) >= 1 && dp_cov == 0)
                reason = "cov2nocov", breakit = true;
            else if (getk_(&mas_dp_) >= 1 && getmax_(&mas_dp_) == 0 && dp_cov >= 1)
                reason = "nocov2cov", breakit = true;
            if (breakit) {
                BASE_t refbase = sequences.get_ref_base(last_block);
                tmppb << region << "\t" << last_block << "\t.\t" << refbase << "\t.\t.\t" << reason
                      << "\tEND=" << (last_block == i || rposit_first - 1 == i ? i : i - 1) << ";BLOCKAVG_"
                      << opts::block_label;
                if (last_block == i) {
                    tmppb << ";PX=" << coverages.get_p(i, i_quals_) << "," << coverages.get_p(i, i_quals_) << ","
                          << coverages.get_p(i, i_quals_) << ",0"
                          << ";VRX=" << (short)coverages.get_force_call(i_vr_, i) << ","
                          << (short)coverages.get_force_call(i_vr_, i) << ","
                          << (short)coverages.get_force_call(i_vr_, i) << ",0"
                          << ";RRX=" << (short)coverages.get_coverage(i_rr_, i) << ","
                          << (short)coverages.get_coverage(i_rr_, i) << "," << (short)coverages.get_coverage(i_rr_, i)
                          << ",0"
                          << ";DPX=" << dp_cov << "," << dp_cov << "," << dp_cov << ",0";
                } else {
                    if (rposit_first - 1 == i) {
                        add_val(&mas_q_, coverages.get_p(i, i_quals_));
                        add_val(&mas_vr_, (short)coverages.get_force_call(i_vr_, i));
                        add_val(&mas_rr_, (short)coverages.get_coverage(i_rr_, i));
                        add_val(&mas_dp_, dp_cov);
                    }
                    tmppb << ";PX=" << getmin_(&mas_q_) << "," << getmax_(&mas_q_) << "," << getmean_(&mas_q_) << ","
                          << getsd_(&mas_q_) << ";VRX=" << getmin_(&mas_vr_) << "," << getmax_(&mas_vr_) << ","
                          << getmean_(&mas_vr_) << "," << getsd_(&mas_vr_) << ";RRX=" << getmin_(&mas_rr_) << ","
                          << getmax_(&mas_rr_) << "," << getmean_(&mas_rr_) << "," << getsd_(&mas_rr_)
                          << ";DPX=" << getmin_(&mas_dp_) << "," << getmax_(&mas_dp_) << "," << getmean_(&mas_dp_)
                          << "," << getsd_(&mas_dp_);
                }
                tmppb << "\tGT:VR:RR:DP:GQ\t0/0:" << getmin_(&mas_vr_) << ":" << getmin_(&mas_rr_) << ":"
                      << getmin_(&mas_dp_) << ":.\n";
                last_block = i;
                reset_blocker(&mas_q_);
                reset_blocker(&mas_vr_);
                reset_blocker(&mas_rr_);
                reset_blocker(&mas_dp_);
            }
            add_val(&mas_q_, coverages.get_p(i, i_quals_));
            add_val(&mas_vr_, (short)coverages.get_force_call(i_vr_, i));
            add_val(&mas_rr_, (short)coverages.get_coverage(i_rr_, i));
            add_val(&mas_dp_, dp_cov);
        }
        output += tmppb.str();
    }
}
inline GENOTYPE indel_genotype(INDEL_EVENT*& p_indel, COV_COUNTER const& coverages) {
    COVERAGE_t n_tmp = 0;
    GENOTYPE genotype;
    memset(&genotype, 0, GENOTYPE_LENGTH * sizeof(char));
    COVERAGE_t total_depth = coverages.get_coverage(COV_COUNTER::AID_DEP, p_indel->var_start_);
    COVERAGE_t total_ref_depth = coverages.get_coverage(COV_COUNTER::AID_REF_DEPTH, p_indel->var_start_);
    if (p_indel->read_count_ == 0) {
        if (total_depth < opts::INDEL_DEPTH_CUTOFF) {
            memcpy(&genotype, "./.", GENOTYPE_LENGTH);
            return genotype;
        }
        n_tmp = total_depth;
    } else {
        if (total_depth < opts::INDEL_DEPTH_CUTOFF) {
            memcpy((char*)&genotype, "1/.", GENOTYPE_LENGTH);
            return genotype;
        }
        n_tmp = p_indel->read_count_ + total_ref_depth;
    }
    if (((double)p_indel->read_count_ / n_tmp) >= opts::INDEL_HOM_VAR_CUTOFF)
        memcpy((char*)&genotype, "1/1", GENOTYPE_LENGTH);
    else if (((double)p_indel->read_count_ / n_tmp) >= opts::INDEL_HET_CUTOFF)
        memcpy((char*)&genotype, "1/0", GENOTYPE_LENGTH);
    else if ((double)total_ref_depth / n_tmp >= opts::INDEL_HOM_VAR_CUTOFF)
        memcpy((char*)&genotype, "0/0", GENOTYPE_LENGTH);
    else if ((double)total_ref_depth / n_tmp >= opts::INDEL_HET_CUTOFF)
        memcpy((char*)&genotype, "0/.", GENOTYPE_LENGTH);
    else
        memcpy((char*)&genotype, "./.", GENOTYPE_LENGTH);
    return genotype;
}
#ifdef STUFF
std::vector<VCF::VAR> snp_evaluate(std::string const& region, REF_ID_t, REF_POS_t pos_,
                                   NASTY_GLOBAL_SNP_HOLDER& events_g, std::ofstream& osnpv, std::ofstream& osnps,
                                   SCF_SEQ const& sequences, COV_COUNTER const& coverages) {
#else
inline void snp_evaluate(std::string const& region, REF_ID_t, REF_POS_t pos_, NASTY_GLOBAL_SNP_HOLDER& events_g,
                         std::ofstream& osnpv, SCF_SEQ const& sequences, COV_COUNTER const& coverages, LIST const& segs,
                         uint32_t piece_i) {
#endif
#ifdef STUFF
    std::vector<VCF::VAR> combined_snp_indel_buffer;
#endif
    std::string site_info = "";
    std::vector<REF_POS_t> wipe_but_should_just_mark;
    for (NASTY_GLOBAL_SNP_HOLDER::SNPS::const_iterator rposit = events_g.snps_.begin(); rposit != events_g.snps_.end();
         rposit++) {
        if (rposit->first < pos_)
            wipe_but_should_just_mark.push_back(rposit->first);
        if (rposit->first < pos_) {
            site_info = ".";
            string read_info_output = "";
            BIG_QUAL_t added_qual = 0;
            BASE_t refbase = sequences.get_ref_base(rposit->first);
            COVERAGE_t refbase_cov = coverages.get_coverage(COV_COUNTER::AS_COV, rposit->first);
            char refenv[READENV];
            memset(refenv, 0, READENV * sizeof(BASE_t));
            if (rposit->first < SREADENV + 1)
                for (int i = 0; i < SREADENV + 1; ++i)
                    refenv[i] = sequences.get_ref_base(rposit->first + i);
            else
                for (int i = 0; i < READENV - 1; ++i)
                    refenv[i] = sequences.get_ref_base(rposit->first + i - (SREADENV));
            uint32_t alternative_reads = events_g.snps_.find(rposit->first)->second.size();
            uint32_t total_coverage =
                refbase_cov + alternative_reads + coverages.get_coverage(COV_COUNTER::AS_DEL, rposit->first);
            if (opts::allnonref || (opts::gvcftest && opts::VRCUTOFF && alternative_reads >= opts::VRCUTOFF &&
                                    (double)(alternative_reads / total_coverage) >= opts::max_alt_frac)) {
            } else if (opts::SLX_SNP_HIGH_MAP_QUAL_MAX_COVERAGE &&
                       total_coverage > opts::SLX_SNP_HIGH_MAP_QUAL_MAX_COVERAGE)
                continue;
            BASE_t high_base = '.';
            COVERAGE_t high_base_cov = 0;
            std::vector<SNP_EVENT> const& tired = events_g.snps_.find(rposit->first)->second;
            std::map<BASE_t, uint32_t> counts;
            std::map<BASE_t, BIG_QUAL_t> quals;
#ifdef STUFF
            SNP_INFOS readgroup_readid_info;
#endif
            for (std::vector<SNP_EVENT>::const_iterator snpit = tired.begin(); snpit != tired.end(); snpit++) {
                counts[snpit->allele_base_]++;
                quals[snpit->allele_base_] += snpit->qual_;
#ifdef STUFF
                readgroup_readid_info[snpit->allele_base_].readids_.push_back(snpit->readid_);
                readgroup_readid_info[snpit->allele_base_].readgroups_[snpit->rg_]++;
#endif
            }
            for (std::map<BASE_t, uint32_t>::iterator bit = counts.begin(); bit != counts.end(); bit++)
                if (bit->second > high_base_cov) {
                    high_base = bit->first;
                    high_base_cov = bit->second;
                }
            std::string filter = "";
            if (high_base_cov > 1) {
                BIG_QUAL_t high_qual = quals[high_base];
                unsigned int equal_majority = 0;
                for (std::map<BASE_t, uint32_t>::iterator bit = counts.begin(); bit != counts.end(); bit++)
                    if (bit->second == high_base_cov)
                        ++equal_majority;
                if (equal_majority > 1) {
                    filter += "equal_majority;";
                    if (!opts::silly_scavenge) {
                        for (std::map<BASE_t, uint32_t>::iterator bit = counts.begin(); bit != counts.end(); bit++)
                            if (bit->second == high_base_cov && quals[bit->first] > high_qual) {
                                high_base = bit->first;
                                high_base_cov = bit->second;
                            }
                    } else {
                        std::map<BASE_t, double> prs;
                        for (std::vector<SNP_EVENT>::const_iterator snpit = tired.begin(); snpit != tired.end();
                             snpit++) {
                            double pr = 1.0 - ill::snp_logit(snpit->qual_, snpit->nqs_, snpit->type_,
                                                             (double)snpit->dist3_ / snpit->read_len_);
                            if (prs.count(snpit->allele_base_) == 0)
                                prs[snpit->allele_base_] = pr;
                            else
                                prs[snpit->allele_base_] *= pr;
                        }
                        for (std::map<BASE_t, double>::iterator bit = prs.begin(); bit != prs.end(); bit++) {
                            uint8_t bin = floor((1.0 - bit->second) * 10);
                            bit->second = (1.0 / (1 + (ill::SLX_err_prior_arr[bin] * ill::prior_err_c /
                                                       ill::SLX_snp_prior_arr[bin] * ill::prior_snp_c)));
                        }
                        double high_pr = prs[high_base];
                        for (std::map<BASE_t, uint32_t>::iterator bit = counts.begin(); bit != counts.end(); bit++)
                            if (bit->second == high_base_cov && prs[bit->first] > high_pr) {
                                high_base = bit->first;
                                high_base_cov = bit->second;
                            }
                    }
                }
            }
            COVERAGE_t var_cov = 0, pos_strand = 0;
            double pr_err_j_prod_sum = 1;
            stringstream rio;
            rio.precision(3);
            for (std::vector<SNP_EVENT>::const_iterator snpit = tired.begin(); snpit != tired.end(); snpit++) {
                double rel_pos = (double)snpit->dist3_ / snpit->read_len_;
                if (snpit->allele_base_ == high_base) {
                    if (snpit->read_strand_)
                        ++pos_strand;
                    added_qual += snpit->qual_;
                    double pr_snp_i_read = ill::snp_logit(snpit->qual_, snpit->nqs_, snpit->type_, rel_pos);
                    double pr_err_i_read = 1 - pr_snp_i_read;
                    pr_err_j_prod_sum *= pr_err_i_read;
                    rio << "(";
                    if (pr_err_i_read > 0.9995)
                        rio << "1.0";
                    else
                        rio << pr_err_i_read;
                    rio << ");";
                    ++var_cov;
                }
            }
            coverages.set_force_call(COV_COUNTER::AS_VR, rposit->first, var_cov);
            coverages.set_force_call(COV_COUNTER::AS_AR, rposit->first, alternative_reads);
            double pr_snp_j_prod_sum = 1 - pr_err_j_prod_sum;
            uint8_t bin = floor(pr_snp_j_prod_sum * 10);
            bin = bin < 10 ? bin : 9;
            double pr_S_j_ERR_c = ill::SLX_err_prior_arr[bin];
            double pr_S_j_SNP_c = ill::SLX_snp_prior_arr[bin];
            double posterior_err = pr_S_j_ERR_c * ill::prior_err_c;
            double posterior_snp = pr_S_j_SNP_c * ill::prior_snp_c;
            double pr_SNP_S_j_c_j = 1.0 / (1 + (posterior_err / posterior_snp));
            if (opts::allnonref || (opts::gvcftest && opts::VRCUTOFF && alternative_reads >= opts::VRCUTOFF &&
                                    (double)(alternative_reads / total_coverage) >= opts::max_alt_frac)) {
            } else if (pr_SNP_S_j_c_j < opts::SNP_MIN_PR)
                continue;
            stringstream rediculous1, rediculous2;
            rediculous1 << pr_S_j_ERR_c;
            string rediculous_s = rediculous1.str();
            rediculous2 << pr_S_j_SNP_c;
            rediculous_s = rediculous2.str();
            COVERAGE_t n_tmp = var_cov + refbase_cov;
#ifdef STUFF
            stringstream ittmp;
            string info_thing = "RG:";
            for (std::map<RG_t, COVERAGE_t>::iterator xit = readgroup_readid_info[high_base].readgroups_.begin();
                 xit != readgroup_readid_info[high_base].readgroups_.end(); xit++)
                ittmp << "\"" << xit->first << "\""
                      << ":" << xit->second << ",";
            ittmp << ";RI:";
            for (std::vector<uint64_t>::iterator xit = readgroup_readid_info[high_base].readids_.begin();
                 xit != readgroup_readid_info[high_base].readids_.end(); xit++)
                ittmp << *xit << ",";
            info_thing += ittmp.str();
#endif
            if (pr_SNP_S_j_c_j < ill::SLX_SNP_CUTOFF)
                filter += "low_snpqual;";
            if (n_tmp > ill::SLX_SNP_HIGH_MAP_QUAL_MAX_COVERAGE)
                filter += "high_coverage;";
            if (n_tmp < ill::SLX_SNP_MIN_COVERAGE)
                filter += "low_coverage;";
            COVERAGE_t neg_strand2 = var_cov - pos_strand;
            double silly1 = (double)neg_strand2 / var_cov;
            double silly2 = (double)pos_strand / var_cov;
            if (n_tmp >= opts::SNP_STRAND_TEST_COV_CUTOFF &&
                ((silly1 < silly2 ? silly1 : silly2) < opts::SNP_STRAND_RAIO_CUTOFF))
                filter += "single_strand;";
            if (var_cov > 0 && var_cov <= opts::SNP_MIN_COV)
                filter += "low_VariantReads;";
            QUAL_t snpq = round(-10 * log10(1 - pr_SNP_S_j_c_j + 0.000001));
            GENOTYPE genotype;
            memset(&genotype, 0, GENOTYPE_LENGTH * sizeof(char));
            if (n_tmp == 0)
                filter += "No_data;";
            else {
                if (((double)var_cov / n_tmp) <= opts::SNP_HET_MIN) {
                    strcpy((char*)&genotype, "0/0");
                    filter += "low_VariantRatio";
                } else if (((double)var_cov / n_tmp) > opts::SNP_HET_MIN &&
                           ((double)var_cov / n_tmp) < opts::SNP_HET_MAX)
                    strcpy((char*)&genotype, "0/1");
                else
                    strcpy((char*)&genotype, "1/1");
            }
            if (n_tmp > 0 && var_cov == 0)
                filter = "No_var";
            if (filter.empty())
                filter = "PASS";
            else if (filter == "equal_majority;")
                filter = "PASS;equal_majority";
            else if (filter[filter.length() - 1] == ';')
                filter = filter.substr(0, filter.length() - 1);
            S_SEQ_t r, a;
            r += refbase;
            a += high_base;
            if (!opts::capturebed ||
                (opts::capturebed && (rposit->first >= segs[piece_i].first && rposit->first <= segs[piece_i].second))) {
                std::string osnpvstr;
                if (opts::gvcftest && last_call_snp != (REF_POS_t)-1) {
                    osnpv.precision(3);
                    output_gvcf(rposit->first, last_call_snp, segs, piece_i, coverages, sequences, region, osnpvstr,
                                COV_COUNTER::SNVP, COV_COUNTER::AS_VR, COV_COUNTER::AS_COV, COV_COUNTER::DUMMY);
                    if (last_call_snp != (REF_POS_t)-1)
                        osnpv << osnpvstr;
                } else if (opts::gvcftestdebug)
                    output_novar(rposit->first, last_call_snp + 1, segs, piece_i, coverages, sequences, region, osnpv,
                                 COV_COUNTER::SNVP, COV_COUNTER::AS_VR, COV_COUNTER::AS_COV, COV_COUNTER::DUMMY);
                osnpv << region << "\t" << rposit->first << "\t.\t" << refbase << "\t" << high_base << "\t"
                      << (short)snpq << "\t" << filter << "\tP=" << pr_SNP_S_j_c_j << "\tGT:VR:RR:DP:GQ\t"
                      << (char*)&genotype << ":" << var_cov << ":" << refbase_cov << ":" << n_tmp << ":.\n";
                if (opts::gvcftest && last_call_snp == (REF_POS_t)-1)
                    osnpv << osnpvstr;
#ifdef STUFF
                combined_snp_indel_buffer.push_back(VCF::VAR(region, rposit->first, r, a, filter, snpq,
                                                             (char*)&genotype, var_cov, alternative_reads, refbase_cov,
                                                             n_tmp, pr_SNP_S_j_c_j, true, info_thing));
#endif
                last_call_snp = rposit->first;
            } else {
#ifdef STUFF
                combined_snp_indel_buffer.push_back(VCF::VAR(region, rposit->first, r, a, filter, snpq,
                                                             (char*)&genotype, var_cov, alternative_reads, refbase_cov,
                                                             n_tmp, pr_SNP_S_j_c_j, false, info_thing));
#endif
                coverages.set_snp_p(rposit->first, pr_SNP_S_j_c_j);
            }
        }
    }
    for (size_t i = 0; i < wipe_but_should_just_mark.size(); i++)
        events_g.snps_.erase(wipe_but_should_just_mark[i]);
#ifdef STUFF
    return combined_snp_indel_buffer;
#endif
}
#ifdef STUFF
inline void output_indel(bool hrc, std::string const& region, INDEL_EVENT*& p_indel, std::ofstream& oindelv,
                         SCF_SEQ const& sequences, COV_COUNTER const& coverages, std::ofstream& oindeli,
                         std::vector<VCF::VAR>& list) {
#else
inline void output_indel(std::string const& region, INDEL_EVENT*& p_indel, std::ofstream& oindelv,
                         SCF_SEQ const& sequences, COV_COUNTER const& coverages, LIST const& segs, uint32_t piece_i) {
#endif
    string filter;
    coverages.set_force_call(COV_COUNTER::AID_VR, p_indel->var_start_, p_indel->read_count_);
    if (p_indel->read_count_ < opts::MIN_VAR_READS)
        filter += "low_VariantReads;";
    COVERAGE_t total_depth = coverages.get_coverage(COV_COUNTER::AID_DEP, p_indel->var_start_);
    if (total_depth < opts::MIN_DEPTH_COVERAGE)
        filter += "low_coverage;";
    if (total_depth > ill::SLX_SNP_HIGH_MAP_QUAL_MAX_COVERAGE)
        filter += "high_coverage;";
    double var_ratio = (double)p_indel->read_count_ / total_depth;
    if (var_ratio < opts::MIN_VAR_RATIO)
        filter += "low_VariantRatio;";
    if (opts::STRAND_DIR_FILTER && !p_indel->simple_strand_test())
        filter += "single_strand;";
    if (p_indel->near_read_end_ratio() > opts::MAX_NEAR_READ_END_RATIO)
        filter += "read_end_ratio;";
    double mpr_strand_dir = p_indel->simple_strand_test() ? 1.0 : 0.0;
    double mean_avnqs = p_indel->mean_avnqs();
    double mean_mapq = p_indel->mean_mapq();
    double mean_var_rate = p_indel->mean_var_rate();
    double mpr_entropy = p_indel->simple_local_entropy(sequences);
    double PR_INDEL_j = 0;
    if (opts::indeltest)
        PR_INDEL_j = ill::indel_logit_wgs(mpr_entropy, mpr_strand_dir, mean_avnqs, mean_mapq, mean_var_rate,
                                          p_indel->read_count_, total_depth);
    else
        PR_INDEL_j = ill::indel_logit(mpr_entropy, mpr_strand_dir, mean_avnqs, p_indel->read_count_, total_depth);
#ifdef STUFF
    stringstream ittmp;
    ittmp << "P_j=" << PR_INDEL_j << ";strand_j=" << mpr_strand_dir << ";avnqs_j=" << mean_avnqs
          << ";simple_entropy=" << mpr_entropy << ";RG:";
    for (std::map<RG_t, COVERAGE_t>::iterator xit = p_indel->readgroups_.begin(); xit != p_indel->readgroups_.end();
         xit++)
        ittmp << "\"" << xit->first << "\""
              << ":" << xit->second << ",";
    ittmp << ";RI:";
    for (std::vector<uint64_t>::iterator xit = p_indel->readids_.begin(); xit != p_indel->readids_.end(); xit++)
        ittmp << *xit << ",";
    string info_thing = ittmp.str();
    info_thing += ";" + p_indel->tmpinfo_;
#endif
    QUAL_t qual = 0;
    double p = 0;
    if (1 - PR_INDEL_j < 0.0000000001)
        qual = 100;
    else
        qual = (double)round(-10.0 * log10(1.0 - PR_INDEL_j));
    p = (double)round(PR_INDEL_j * 10000.0) / 10000.0;
    if (opts::allnonref ||
        (opts::gvcftest && opts::VRCUTOFF && p_indel->read_count_ >= opts::VRCUTOFF &&
         (double)p_indel->read_count_ / coverages.get_coverage(COV_COUNTER::AID_DEP, p_indel->var_start_) >=
             opts::max_alt_frac)) {
    } else if (qual == 0 || p <= opts::INDEL_PR_CUTOFF ||
               (opts::SLX_SNP_HIGH_MAP_QUAL_MAX_COVERAGE &&
                coverages.get_coverage(COV_COUNTER::AID_DEP, p_indel->var_start_) >
                    opts::SLX_SNP_HIGH_MAP_QUAL_MAX_COVERAGE))
        return;
    GENOTYPE genotype = indel_genotype(p_indel, coverages);
    string alt, ref;
    if (p_indel->isdel_) {
        alt = sequences.get_ref_base(p_indel->var_start_);
        ref = sequences.get_ref_seq_range(p_indel->var_start_, p_indel->var_len_ + 1);
    } else {
        ref = sequences.get_ref_base(p_indel->var_start_);
        alt = ref;
        alt += (strchr(p_indel->indel_.data(), 'I') + 1);
    }
    if (p < 0.5)
        filter += "low_qual;";
    if (filter.empty())
        filter += "PASS";
    else if (filter[filter.length() - 1] == ';')
        filter = filter.substr(0, filter.length() - 1);
    if (!opts::capturebed || (opts::capturebed && (p_indel->var_start_ >= segs[piece_i].first &&
                                                   p_indel->var_start_ <= segs[piece_i].second))) {
        std::string oindelvstr;
        if (opts::gvcftest && last_call_indel != (REF_POS_t)-1) {
            oindelv.precision(3);
            output_gvcf(p_indel->var_start_, last_call_indel, segs, piece_i, coverages, sequences, region, oindelvstr,
                        COV_COUNTER::INDELP, COV_COUNTER::AID_VR, COV_COUNTER::AID_REF_DEPTH, COV_COUNTER::AID_DEP);
            if (last_call_indel != (REF_POS_t)-1)
                oindelv << oindelvstr;
        } else if (opts::gvcftestdebug)
            output_novar(p_indel->var_start_, last_call_indel + 1, segs, piece_i, coverages, sequences, region, oindelv,
                         COV_COUNTER::INDELP, COV_COUNTER::AID_VR, COV_COUNTER::AID_REF_DEPTH, COV_COUNTER::AID_DEP);
        oindelv << region << "\t" << p_indel->var_start_ << "\t.\t" << ref << "\t" << alt << "\t" << (short)qual << "\t"
                << filter << "\tP=";
        if (opts::compatibility && p == 1.0)
            oindelv << "1.0";
        else
            oindelv << p;
        oindelv << "\tGT:VR:RR:DP:GQ\t" << (char*)&genotype << ":" << p_indel->read_count_ << ":"
                << coverages.get_coverage(COV_COUNTER::AID_REF_DEPTH, p_indel->var_start_) << ":"
                << coverages.get_coverage(COV_COUNTER::AID_DEP, p_indel->var_start_) << ":.\n";
        if (opts::gvcftest && last_call_indel == (REF_POS_t)-1)
            oindelv << oindelvstr;
#ifdef STUFF
        list.push_back(
            VCF::VAR(region, p_indel->var_start_, ref, alt, filter, qual, (char*)&genotype, p_indel->read_count_,
                     coverages.get_force_call(COV_COUNTER::AID_AR, p_indel->var_start_),
                     coverages.get_coverage(COV_COUNTER::AID_REF_DEPTH, p_indel->var_start_),
                     coverages.get_coverage(COV_COUNTER::AID_DEP, p_indel->var_start_), p, true, info_thing));
#endif
        last_call_indel = p_indel->var_start_;
    } else {
#ifdef STUFF
        list.push_back(
            VCF::VAR(region, p_indel->var_start_, ref, alt, filter, qual, (char*)&genotype, p_indel->read_count_,
                     coverages.get_force_call(COV_COUNTER::AID_AR, p_indel->var_start_),
                     coverages.get_coverage(COV_COUNTER::AID_REF_DEPTH, p_indel->var_start_),
                     coverages.get_coverage(COV_COUNTER::AID_DEP, p_indel->var_start_), p, false, info_thing));
#endif
        coverages.set_indel_p(p_indel->var_start_, p);
    }
}
#ifdef STUFF
std::vector<VCF::VAR> print_indel_buffer(std::string const& region, REF_POS_t bpos, NASTY_GLOBAL_SNP_HOLDER& events_g,
                                         std::ofstream& oindelv, SCF_SEQ const& sequences, COV_COUNTER const& coverages,
                                         std::ofstream& oindeli) {
    std::vector<VCF::VAR> combined_snp_indel_buffer;
#else
void print_indel_buffer(std::string const& region, REF_POS_t bpos, NASTY_GLOBAL_SNP_HOLDER& events_g,
                        std::ofstream& oindelv, SCF_SEQ const& sequences, COV_COUNTER const& coverages,
                        LIST const& segs, uint32_t piece_i) {
#endif
    NASTY_GLOBAL_SNP_HOLDER::INDELS& it = events_g.indel_;
    std::vector<REF_POS_t> keys_to_wipe;
    INDEL_EVENT* previous_indel = 0;
    for (NASTY_GLOBAL_SNP_HOLDER::INDELS::iterator xit2 = it.begin(); xit2 != it.end(); xit2++) {
        for (NASTY_GLOBAL_SNP_HOLDER::INDEL_GROUP::iterator xit = xit2->second.begin(); xit != xit2->second.end();
             xit++) {
            if (previous_indel) {
                if (it[xit2->first][xit->first].var_start_ == previous_indel->var_start_) {
                    if (previous_indel->read_count_ < it[xit2->first][xit->first].read_count_) {
#ifdef STUFF
                        output_indel(false, region, previous_indel, oindelv, sequences, coverages, oindeli,
                                     combined_snp_indel_buffer);
#else
                        if (opts::allindels)
                            output_indel(region, previous_indel, oindelv, sequences, coverages, segs, piece_i);
#endif
                        delete previous_indel;
                        previous_indel = new INDEL_EVENT(it[xit2->first][xit->first]);
                        coverages.add_force_call(COV_COUNTER::AID_AR, it[xit2->first][xit->first].var_start_);
                    }
                    keys_to_wipe.push_back(xit2->first);
                    continue;
                } else {
#ifdef STUFF
                    output_indel(true, region, previous_indel, oindelv, sequences, coverages, oindeli,
                                 combined_snp_indel_buffer);
#else
                    output_indel(region, previous_indel, oindelv, sequences, coverages, segs, piece_i);
#endif
                    delete previous_indel;
                    previous_indel = 0;
                }
            }
            cout.flush();
            if (it[xit2->first][xit->first].var_start_ < bpos) {
                keys_to_wipe.push_back(xit2->first);
                if (previous_indel)
                    delete previous_indel;
                previous_indel = new INDEL_EVENT(it[xit2->first][xit->first]);
            }
        }
    }
    if (previous_indel != 0) {
#ifdef STUFF
        output_indel(true, region, previous_indel, oindelv, sequences, coverages, oindeli, combined_snp_indel_buffer);
#else
        output_indel(region, previous_indel, oindelv, sequences, coverages, segs, piece_i);
#endif
        delete previous_indel;
    }
    for (std::vector<REF_POS_t>::iterator wit = keys_to_wipe.begin(); wit != keys_to_wipe.end(); wit++)
        events_g.indel_.erase(*wit);
#ifdef STUFF
    return combined_snp_indel_buffer;
#endif
}
}
bool check_file_okay(char const* f) {
    fstream xx;
    xx.open(f, std::ios::in | std::ios::binary);
    if (!xx) {
        cerr << "problem opening file\n";
        return false;
    }
    xx.seekg(0, ios::end);
    unsigned long long size = xx.tellg();
    char* memblock = new char[size];
    xx.seekg(0, ios::beg);
    xx.read(memblock, size);
    xx.close();
    int probs = 0;
    for (uint64_t i = 0; i < size; ++i) {
        if (memblock[i] == '\0')
            ++probs;
    }
    delete[] memblock;
    return probs == 0;
}
inline bool nastyness(bool a, char const* const pfx, char const* const region) {
    std::string checkpoint(pfx);
    checkpoint + region;
    checkpoint += "_seq";
    checkpoint += string(region);
    checkpoint += ".done";
    if (a) {
        if (utils::isregfile(checkpoint)) {
            if (!opts::force) {
                cerr << "region " << region << " has already been processed, will not reporcess without force ("
                     << checkpoint << ")\n";
                return false;
            } else {
                utils::removefile(checkpoint);
                cerr << "region " << region << " will be re-run (" << checkpoint << ")\n";
            }
        } else
            cerr << "no checkpoint file present\n";
    } else {
        if (opts::checkpoint) {
            cerr << "generating checkpoint file\n";
            utils::touchfile(checkpoint);
        } else
            cerr << "not generating checkpoint file\n";
    }
    return true;
}
void do_it(LIST& segs, char const* region, char const* bf, const char* refseq, QUAL_t min_map_qual,
           QUAL_t min_map_qual2, char const* pfx) {
    if (!nastyness(true, pfx, region))
        return;
    string tmps(region);
    uint32_t piece_i = 0;
    std::ofstream oindelv;
    oindelv.open((std::string(pfx) + "_seq" + tmps + "_indel.vcf").data(), std::fstream::out | std::fstream::trunc);
    std::ofstream osnpv;
    osnpv.open((std::string(pfx) + "_seq" + tmps + "_snp.vcf").data(), std::fstream::out | std::fstream::trunc);
    KNOWN_INDELS ki;
    if (opts::known_indels)
        ki = get_known(tmps);
#ifdef STUFF
    std::ofstream ov;
    ov.open((std::string(pfx) + "_seq" + tmps + ".vcf").data(), std::fstream::out | std::fstream::trunc);
#endif
    if (!oindelv)
        cerr << "unable to open output indel vcf file\n", exit(1);
    if (!osnpv)
        cerr << "unable to open output snp vcf file\n", exit(1);
#ifdef STUFF
    if (!ov)
        cerr << "unable to open output vcf file\n", exit(1);
#endif
    std::ofstream osnps, oindeli;
    if (opts::snp_file) {
        osnps.open((std::string(pfx) + "_seq" + tmps + "_snp.snp").data());
        osnps.precision(3);
        if (!osnps)
            cerr << "unable to open output snp snp file\n", exit(1);
        oindeli.open((std::string(pfx) + "_seq" + tmps + "_indel.indel").data());
        oindeli.precision(3);
        if (!oindeli)
            cerr << "unable to open output indel indel file\n", exit(1);
    }
    NASTY_GLOBAL_SNP_HOLDER events_g;
    uint32_t dodgy_reads = 0, missing_nm_tag = 0;
    SCF_SEQ sequences(refseq, region);
    COV_COUNTER coverages(sequences.get_length());
    std::cerr << "fetched sequence for " << region << " of length " << sequences.get_length() << "\n";
    STUFF bam(bf);
    std::vector<std::string> piece;
    if (segs.size() == 0 && !opts::capturebed) {
        piece.push_back(region);
        segs.push_back(make_pair(1, sequences.get_length()));
    } else {
        for (unsigned int i = 0; i < segs.size(); ++i) {
            std::stringstream x;
            x << region << ":" << segs[i].first << "-" << segs[i].second;
            piece.push_back(x.str());
        }
    }
    cerr << "fetching reads for " << region << " (" << segs.size() << ")\n";
    bam1_t* b = bam.alignment;
    uint64_t read_id = 0;
    timespec nsstart, nsend;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &nsstart);
    int ret;
    std::vector<uint8_t> cigars_per_string;
    uint8_t cigars = 0;
    for (piece_i = 0; piece_i < piece.size(); ++piece_i) {
        last_call_snp = segs[piece_i].first - 1;
        last_call_indel = segs[piece_i].first - 1;
        string const& what = piece[piece_i];
        hts_itr_t* iter = bam.get_iter(what.data());
        {
            while ((ret = sam_itr_next(bam.in, iter, b)) >= 0) {
                ++read_id;
                if (opts::rate) {
                    if (read_id % 10000 == 0) {
                        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &nsend);
                        long nsd = nsend.tv_nsec - nsstart.tv_nsec;
                        cerr << "processed " << read_id << " reads @ " << b->core.pos << " ("
                             << (double)10000000000000 / abs(nsd) << " reads/s)\r";
                        nsstart = nsend;
                    }
                }
                if (cigars) {
                    cigars_per_string.push_back(cigars);
                    cigars = 0;
                }
                uint32_t* cit = bam_get_cigar(b);
                bool strand = bam_is_rev(b) == 0;
                RG_t rg;
                int32_t score_as = -1, nm = -1;
                uint16_t fl = b->core.flag;
                if (fl & BAM_FUNMAP || fl & BAM_FSECONDARY || fl & BAM_FQCFAIL || fl & BAM_FDUP)
                    continue;
                if (coverages.is_toxic(COV_COUNTER::AID_DEP, b->core.pos) &&
                    coverages.is_toxic(COV_COUNTER::AID_DEP, (b->core.pos + (b->core.l_qseq >> 1))) &&
                    coverages.is_toxic(COV_COUNTER::AID_DEP, b->core.pos + b->core.l_qseq)) {
                    continue;
                } else if (b->core.n_cigar == 0)
                    continue;
                else if (b->core.qual < min_map_qual)
                    continue;
                uint8_t* s = bam_get_aux(b);
                while (s < b->data + b->l_data) {
                    uint8_t type, key[3];
                    key[2] = '\0';
                    key[0] = s[0];
                    key[1] = s[1];
                    s += 2;
                    type = *s;
                    ++s;
                    if (memcmp(key, (void*)"NM", 2) == 0) {
                        assert(type != 'C' || type != 'c');
                        nm = *s;
                    } else if (memcmp(key, (void*)"AS", 2) == 0)
                        score_as = *s;
                    else if (memcmp(key, (void*)"RG", 2) == 0)
                        rg = (char*)s;
                    switch (type) {
                    case 'A':
                    case 'C':
                    case 'c':
                        ++s;
                        break;
                    case 'S':
                    case 's':
                        s += 2;
                        break;
                    case 'I':
                    case 'i':
                    case 'f':
                        s += 4;
                        break;
                    case 'd':
                        s += 8;
                        break;
                    case 'Z':
                    case 'H':
                        while (*s)
                            s++;
                        ++s;
                        break;
                    case 'B':
                        break;
                    default:
                        assert(0);
                        break;
                    }
                }
                if (nm == -1) {
                    ++missing_nm_tag;
                    continue;
                }
                A_READS_INDEL_LIST this_reads_indels;
                int32_t aid_offset = 0;
                REF_POS_t aid_pos = b->core.pos + 1;
                REF_POS_t pos_ = aid_pos, as_filtered = 0,
                          leading_sclip = ((*cit & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP)
                                              ? ((*cit & ~(BAM_CIGAR_MASK)) >> BAM_CIGAR_SHIFT)
                                              : 0,
                          tailing_sclip = ((*(cit + b->core.n_cigar - 1) & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP)
                                              ? ((*(cit + b->core.n_cigar - 1) & ~(BAM_CIGAR_MASK)) >> BAM_CIGAR_SHIFT)
                                              : 0;
                uint8_t as_nt_snp = 0, as_nt_ins = 0, as_nt_del = 0, as_nt_sclip = 0;
                bool aid_prev_ins = false, aid_first_nt = true;
                for (size_t i = 0; i < b->core.n_cigar; ++i, ++cit) {
                    ++cigars;
                    int32_t c_op = cigar_op(cit);
                    REF_POS_t c_len = cigar_len(cit);
                    if (b->core.qual >= min_map_qual2) {
                        switch (c_op) {
                        case BAM_CMATCH:
                            break;
                        case BAM_CINS:
                            as_nt_ins += c_len;
                            break;
                        case BAM_CDEL:
                            as_nt_del += c_len;
                            break;
                        case BAM_CSOFT_CLIP:
                            as_nt_sclip += c_len;
                            break;
                        default:
                            assert(0);
                            break;
                        }
                    }
                    switch (c_op) {
                    case BAM_CMATCH:
                    case BAM_CREF_SKIP:
                    case BAM_CEQUAL:
                    case BAM_CDIFF: {
                        for (uint32_t i3 = aid_pos - 1; i3 < aid_pos + c_len - 1; ++i3) {
                            if (aid_first_nt) {
                                aid_first_nt = false;
                                continue;
                            }
                            if (coverages.seen(i3) == false) {
                                coverages.set_coverage(COV_COUNTER::AID_DEP, i3, 1);
                                coverages.set_coverage(COV_COUNTER::AID_REF_DEPTH, i3, 1);
                            } else {
                                coverages.add_coverage(COV_COUNTER::AID_DEP, i3);
                                coverages.add_coverage(COV_COUNTER::AID_REF_DEPTH, i3);
                            }
                            if (aid_prev_ins) {
                                coverages.subtract_coverage(COV_COUNTER::AID_DEP, i3);
                                coverages.subtract_coverage(COV_COUNTER::AID_REF_DEPTH, i3);
                                aid_prev_ins = false;
                            }
                        }
                        aid_offset += c_len, aid_pos += c_len;
                    } break;
                    default: {
                        if (c_op == BAM_CDEL) {
                            for (uint32_t i4 = aid_pos - 1; i4 < aid_pos + c_len - 1; ++i4) {
                                if (coverages.seen(i4) == false) {
                                    coverages.set_coverage(COV_COUNTER::AID_DEP, i4, 1);
                                    coverages.set_coverage(COV_COUNTER::AID_REF_DEPTH, i4, 0);
                                } else {
                                    coverages.add_coverage(COV_COUNTER::AID_DEP, i4);
                                }
                            }
                            this_reads_indels.push_back(A_READS_INDEL(b->core.tid, aid_pos - 1, string(""),
                                                                      leading_sclip, tailing_sclip, aid_offset - 1,
                                                                      c_len));
                            aid_pos += c_len;
                            aid_prev_ins = false;
                        } else if (c_op == BAM_CPAD) {
                            aid_pos += c_len;
                            aid_prev_ins = false;
                        } else if (c_op == BAM_CINS) {
                            if (coverages.seen(aid_pos - 1) == false) {
                                coverages.set_coverage(COV_COUNTER::AID_DEP, aid_pos - 1, 1);
                                coverages.set_coverage(COV_COUNTER::AID_REF_DEPTH, aid_pos - 1, 0);
                            } else {
                                coverages.add_coverage(COV_COUNTER::AID_DEP, aid_pos - 1);
                            }
                            std::string name;
                            for (size_t i5 = aid_offset + leading_sclip; i5 < aid_offset + leading_sclip + c_len; ++i5)
                                name += seq_nt16_str[bam_seqi(bam_get_seq(b), i5)];
                            this_reads_indels.push_back(A_READS_INDEL(b->core.tid, aid_pos - 1, name, leading_sclip,
                                                                      tailing_sclip, aid_offset - 1, c_len));
                            aid_offset += c_len;
                            aid_prev_ins = true;
                            aid_first_nt = false;
                        } else if (c_op == BAM_CSOFT_CLIP || c_op == BAM_CHARD_CLIP) {
                        } else
                            assert(0);
                    } break;
                    }
                }
                if (this_reads_indels.size() > 0) {
                    for (A_READS_INDEL_LIST::iterator it = this_reads_indels.begin(); it != this_reads_indels.end();
                         it++) {
                        uint32_t var_count = 0, var_count2 = 0, read_pos = 0, ref_pos = 0;
                        cit = bam_get_cigar(b);
                        uint16_t cigar_i = 0;
                        if (leading_sclip) {
                            ++cit;
                            cigar_i = 1;
                        }
                        NEARBY_VAR_SITES nearby_var_sites;
                        LENGTH_t indel_length = 0;
                        for (; cigar_i < b->core.n_cigar; ++cigar_i, ++cit) {
                            ++cigars;
                            int32_t c_op = cigar_op(cit);
                            REF_POS_t c_len = cigar_len(cit);
                            if (cigar_i == b->core.n_cigar - 1 && c_op == BAM_CSOFT_CLIP)
                                break;
                            if (c_op == BAM_CINS) {
                                ++var_count;
                                var_count2 += c_len;
                                read_pos += c_len;
                                indel_length += (LENGTH_t)c_len;
                            } else if (c_op == BAM_CDEL) {
                                ++var_count;
                                var_count2 += c_len;
                                ref_pos += c_len;
                                indel_length += (LENGTH_t)c_len;
                            } else if (c_op == BAM_CPAD)
                                ref_pos += c_len;
                            else if (c_op == BAM_CMATCH || c_op == BAM_CEQUAL || c_op == BAM_CDIFF) {
                                for (READ_POS_t j = read_pos; j < read_pos + c_len; ++j) {
                                    REF_POS_t const check = get_read_ref_pos(b) + ref_pos + 1;
                                    if ((c_op != BAM_CEQUAL &&
                                         get_read_base(b, j + leading_sclip) != sequences.get_ref_base(check)) ||
                                        c_op == BAM_CDIFF) {
                                        ++var_count;
                                        ++var_count2;
                                    }
                                    ++ref_pos;
                                }
                                read_pos += c_len;
                            } else if (c_op == BAM_CREF_SKIP) {
                                for (READ_POS_t j = read_pos; j < read_pos + c_len; ++j)
                                    ++ref_pos;
                                read_pos += c_len;
                            } else {
                                cigar2asci(b);
                                assert(0);
                            }
                        }
                        LENGTH_t rlen = b->core.l_qseq - it->lsclip - it->tsclip;
                        REF_POS_t read_end_pos = 0, read_ave_nqs = 0, upper = 0;
                        int32_t lower = 0;
                        if (it->seq_.empty()) {
                            read_end_pos = it->offset + 1;
                            lower = (it->offset < (NEAR_END + 1)) ? 0 : it->offset - (NEAR_END + 1);
                            int32_t const x1 = (it->offset + (NEAR_END + 1)), x2 = rlen - 1;
                            upper = x1 > x2 ? rlen - 1 : x1;
                        } else {
                            read_end_pos = it->offset + it->var_len_ + 1;
                            int32_t const x1 = (it->offset - (NEAR_END + 1));
                            lower = x1 < 0 ? 0 : x1;
                            upper = it->offset + (NEAR_END + 1) + it->var_len_ > rlen - 1
                                        ? rlen - 1
                                        : it->offset + (NEAR_END + 1) + it->var_len_;
                        }
                        uint32_t qualsum = 0;
                        QUAL_t* quals = bam_get_qual(b) + it->lsclip;
                        for (uint32_t i = lower; i < upper + 1; ++i)
                            qualsum += quals[i];
                        read_ave_nqs = (float)qualsum / (upper - lower + 1);
                        bool near_read_end = ((it->offset < (NEAR_END + 1)) || (rlen - read_end_pos < (NEAR_END + 2)));
                        double var_rate_gap_and_mismatch_ = round(((double)var_count / rlen) * 10000) / 10000;
                        events_g.add_indel(it->var_start_,
                                           INDEL_EVENT(it->ref_tid_, it->var_start_, it->var_len_, strand, !(strand),
                                                       near_read_end, read_ave_nqs, b->core.qual, it->silly_key(),
                                                       var_rate_gap_and_mismatch_, it->seq_.empty(), read_id, rg));
                    }
                }
                if (b->core.qual < min_map_qual2)
                    continue;
                uint8_t as_nt_gap = as_nt_ins + as_nt_del, as_nt_sub = nm - as_nt_gap;
                float as_sub_percent = 10000 * ((float)as_nt_sub / (b->core.l_qseq - as_nt_sclip)) / 100.0;
                float as_gap_percent = 10000 * ((float)as_nt_gap / (b->core.l_qseq - as_nt_sclip)) / 100.0;
                if (as_sub_percent > ill::SLX_SNP_MAX_SUB || as_gap_percent > ill::SLX_SNP_MAX_INDEL) {
                    ++as_filtered;
                    continue;
                }
                REF_POS_t as_tplace = pos_, as_qplace = 1, as_span = 0;
                int32_t as_cql = 0;
                cit = bam_get_cigar(b);
                A_READS_SNP_LIST this_reads_snps;
                for (size_t i = 0; i < b->core.n_cigar; ++i, ++cit) {
                    ++cigars;
                    int32_t c_op = (*cit & BAM_CIGAR_MASK);
                    REF_POS_t c_len = ((*cit & ~(BAM_CIGAR_MASK)) >> BAM_CIGAR_SHIFT);
                    if (c_op == BAM_CMATCH) {
                        as_cql += c_len;
                        as_span += c_op;
                        for (uint32_t i6 = 0; i6 < c_len; ++i6) {
                            char refbase = sequences.get_ref_base(as_tplace);
                            char allele = seq_nt16_str[bam_seqi(bam_get_seq(b), as_qplace - 1)];
                            char qbase = ::toupper(allele);
                            if (as_nt_sub == 0 || refbase == qbase)
                                coverages.add_coverage(COV_COUNTER::AS_COV, as_tplace);
                            else {
                                uint8_t qual = *(bam_get_qual(b) + as_qplace - 1);
                                uint32_t as_3p = strand ? b->core.l_qseq - as_qplace : as_qplace - 1;
                                bool nqs = as_nqs(as_qplace - 1, as_3p, bam_get_qual(b), b->core.l_qseq);
                                string as_readenv("");
                                if (as_qplace == 1) {
                                    as_readenv += qbase;
                                    for (size_t i7 = as_qplace; i7 < as_qplace + 6; ++i7)
                                        as_readenv += ::tolower(seq_nt16_str[bam_seqi(bam_get_seq(b), i7)]);
                                } else if (as_qplace < 7) {
                                    for (size_t i7 = 0; i7 < as_qplace - 1; ++i7)
                                        as_readenv += ::tolower(seq_nt16_str[bam_seqi(bam_get_seq(b), i7)]);
                                    as_readenv += qbase;
                                    for (size_t i7 = as_qplace; i7 < as_qplace + 6; ++i7)
                                        as_readenv += ::tolower(seq_nt16_str[bam_seqi(bam_get_seq(b), i7)]);
                                } else {
                                    for (size_t i7 = as_qplace - 7; i7 < as_qplace - 1; ++i7)
                                        as_readenv += ::tolower(seq_nt16_str[bam_seqi(bam_get_seq(b), i7)]);
                                    as_readenv += qbase;
                                    for (size_t i7 = as_qplace;
                                         i7 < (b->core.l_qseq - as_qplace > 6 ? as_qplace + 6 : b->core.l_qseq); ++i7)
                                        as_readenv += ::tolower(seq_nt16_str[bam_seqi(bam_get_seq(b), i7)]);
                                }
                                if (!strand) {
                                    for (uint32_t i7 = 0; i7 < (as_readenv.length() >> 1); ++i7) {
                                        char tmp = as_readenv[as_readenv.length() - i7 - 1];
                                        cuidado::switch_it(as_readenv[i7], as_readenv, as_readenv.length() - i7 - 1);
                                        cuidado::switch_it(tmp, as_readenv, i7);
                                    }
                                    if (as_readenv.length() >= 3 && as_readenv.length() % 2)
                                        cuidado::switch_it(as_readenv[(as_readenv.length() >> 1)], as_readenv,
                                                           (as_readenv.length() >> 1));
                                }
                                if (as_readenv.empty())
                                    as_readenv = "";
                                this_reads_snps.insert(make_pair(
                                    as_tplace, SNP_EVENT(allele, refbase, qual, as_qplace, as_3p, nqs,
                                                         as_readenv.data(), bam_get_qname(b), strand, b->core.l_qseq,
                                                         score_as, as_sub_percent, as_gap_percent, read_id, rg)));
                                ++as_nt_snp;
                            }
                            as_tplace += 1, as_qplace += 1;
                        }
                    } else if (c_op == BAM_CINS || c_op == BAM_CSOFT_CLIP)
                        as_qplace += c_len, as_cql += c_len;
                    else if (c_op == BAM_CDEL) {
                        for (POS_t i7 = 0; i7 < c_len; ++i7) {
                            coverages.add_coverage(COV_COUNTER::AS_DEL, as_tplace);
                            ++as_tplace;
                        }
                        as_span += c_len;
                    } else
                        assert(0);
                }
                if (as_nt_sub != as_nt_snp) {
                    ++dodgy_reads;
                    continue;
                }
                assert(as_cql == b->core.l_qseq);
                for (A_READS_SNP_LIST::iterator it = this_reads_snps.begin(); it != this_reads_snps.end(); it++) {
                    char refb = it->second.ref_base_;
                    char snpb = it->second.allele_base_;
                    if (refb == 'N' || snpb == 'N')
                        continue;
                    bool neighbour = false;
                    uint32_t pos = it->first;
                    for (uint32_t i = pos - 2; i < pos + 3; ++i) {
                        if (i == pos)
                            continue;
                        else if (this_reads_snps.count((i)) != 0)
                            neighbour = true;
                    }
                    if (neighbour) {
                        bool snpswap = false;
                        int should_the_order_matter[] = {1, -1, +2, -2};
                        for (uint32_t i = 0; i < sizeof(should_the_order_matter) / sizeof(should_the_order_matter[0]);
                             ++i) {
                            if (this_reads_snps.count(pos + should_the_order_matter[i]) != 0 &&
                                refb == this_reads_snps.find(pos + should_the_order_matter[i])->second.allele_base_ &&
                                snpb == sequences.get_ref_base(pos + should_the_order_matter[i])) {
                                snpswap = true;
                                break;
                            }
                        }
                        if (snpswap)
                            it->second.type_ = SNP_EVENT::SWAP_;
                        else if ((this_reads_snps.count(pos + 1) != 0 &&
                                  this_reads_snps.find(pos + 1)->second.allele_base_ != 'N') ||
                                 (this_reads_snps.count(pos - 1) != 0 &&
                                  this_reads_snps.find(pos - 1)->second.allele_base_ != 'N'))
                            it->second.type_ = SNP_EVENT::MNP_;
                        else
                            it->second.type_ = SNP_EVENT::SNP_;
                    } else
                        it->second.type_ = SNP_EVENT::SNP_;
                    events_g.add_snp(it->first, it->second);
                }
#ifdef STUFF
                std::vector<VCF::VAR> snps = snp_evaluate(region, b->core.tid, pos_, events_g, osnpv, osnps, sequences,
                                                          coverages),
                                      indels;
                if (events_g.indel_.size() > 0)
                    indels = print_indel_buffer(region, pos_ - 1, events_g, oindelv, sequences, coverages, oindeli);
#else
                snp_evaluate(region, b->core.tid, pos_, events_g, osnpv, sequences, coverages, segs, piece_i);
                if (events_g.indel_.size() > 0)
                    print_indel_buffer(region, pos_ - 1, events_g, oindelv, sequences, coverages, segs, piece_i);
#endif
#ifdef STUFF
                if (snps.size() > 0 && indels.size() > 0)
                    ov << "\n------------------------------\n";
                else if (snps.size() > 0 || indels.size() > 0)
                    ov << "\n-----------------\n";
                for (unsigned int i = 0; i < snps.size(); ++i)
                    ov << snps[i];
                for (unsigned int i = 0; i < indels.size(); ++i)
                    ov << indels[i];
#endif
            }
        }
        if (events_g.snps_.size() > 0)
            snp_evaluate(region, b->core.tid, ((REF_POS_t)-1), events_g, osnpv, sequences, coverages, segs, piece_i);
        if (events_g.indel_.size() > 0)
            print_indel_buffer(region, ((REF_POS_t)-1), events_g, oindelv, sequences, coverages, segs, piece_i);
        if (opts::gvcftest) {
            if (last_call_snp != (REF_POS_t)-1) {
                std::string osnpvstr;
                last_call_snp = last_call_snp < segs[piece_i].first ? segs[piece_i].first - 1 : last_call_snp;
                output_gvcf(segs[piece_i].second + 1, last_call_snp, segs, piece_i, coverages, sequences, region,
                            osnpvstr, COV_COUNTER::SNVP, COV_COUNTER::AS_VR, COV_COUNTER::AS_COV, COV_COUNTER::DUMMY);
                osnpv << osnpvstr;
            }
            if (last_call_indel != (REF_POS_t)-1) {
                std::string oindelvstr;
                last_call_indel = last_call_indel < segs[piece_i].first ? segs[piece_i].first - 1 : last_call_indel;
                output_gvcf(segs[piece_i].second + 1, last_call_indel, segs, piece_i, coverages, sequences, region,
                            oindelvstr, COV_COUNTER::INDELP, COV_COUNTER::AID_VR, COV_COUNTER::AID_REF_DEPTH,
                            COV_COUNTER::AID_DEP);
                oindelv << oindelvstr;
            }
        } else if (opts::gvcftestdebug) {
            last_call_snp = last_call_snp < segs[piece_i].first ? segs[piece_i].first : last_call_snp + 1;
            output_novar(segs[piece_i].second + 1, last_call_snp, segs, piece_i, coverages, sequences, region, osnpv,
                         COV_COUNTER::SNVP, COV_COUNTER::AS_VR, COV_COUNTER::AS_COV, COV_COUNTER::DUMMY);
            last_call_indel = last_call_indel < segs[piece_i].first ? segs[piece_i].first : last_call_indel + 1;
            output_novar(segs[piece_i].second + 1, last_call_indel, segs, piece_i, coverages, sequences, region,
                         oindelv, COV_COUNTER::INDELP, COV_COUNTER::AID_VR, COV_COUNTER::AID_REF_DEPTH,
                         COV_COUNTER::AID_DEP);
        }
    }
    assert(events_g.snps_.size() == 0);
    assert(events_g.indel_.size() == 0);
    if (opts::rate)
        cerr << "\n";
    cerr << "ignored " << dodgy_reads << " reads\n";
    oindelv.close();
    osnpv.close();
    osnps.close();
#ifdef STUFF
    ov.close();
#endif
    oindeli.close();
    if (check_file_okay((std::string(pfx) + "_seq" + tmps + "_indel.vcf").data()) &&
        check_file_okay((std::string(pfx) + "_seq" + tmps + "_snp.vcf").data())) {
        nastyness(false, pfx, region);
        return;
    } else
        cerr << "there is a serious consistency issue with output file generation\n", exit(1);
}
static unsigned int restart_pos = 0;
static int restart_int[] = {10, 60, 300, 900, 1800};
void sig_action_function(int, siginfo_t* info, void*) {
    std::cout << "\n---\nMessage received from child : '" << (char*)info->si_value.sival_ptr << "'\n---\n" << std::endl;
    if (strcmp((char*)info->si_value.sival_ptr, "All done") == 0)
        std::cout << "\n---\nParent fork shutting down cos we're done\n", exit(0);
    else
        std::cout << "Parent fork shutting down cos no idea what's going on\n", exit(1);
}
typedef struct _sig_ucontext {
    unsigned long uc_flags;
    struct ucontext* uc_link;
    stack_t uc_stack;
    struct sigcontext uc_mcontext;
    sigset_t uc_sigmask;
} sig_ucontext_t;
void crit_err_hdlr(int sig_num, siginfo_t* info, void* ucontext) {
    void* array[50];
    void* caller_address;
    char** messages;
    int size, i;
    sig_ucontext_t* uc;
    uc = (sig_ucontext_t*)ucontext;
    caller_address = (void*)uc->uc_mcontext.rip;
    fprintf(stderr, "signal %d (%s), address is %p from %p\n", sig_num, strsignal(sig_num), info->si_addr,
            (void*)caller_address);
    size = backtrace(array, 50);
    array[1] = caller_address;
    messages = backtrace_symbols(array, size);
    for (i = 1; i < size && messages != NULL; ++i) {
        fprintf(stderr, "[bt]: (%d) %s\n", i, messages[i]);
    }
    free(messages);
    exit(1);
}
void do_it_all(LISTOREGIONS& listx, char const* a, char const* b, char const* pfx, QUAL_t mmqual, QUAL_t mmqual2) {
#ifdef _OPENMP
    cerr << "using openmp\n";
#pragma omp parallel for ordered schedule(static)
#else
    cerr << "not threaded\n";
#endif
    for (LISTOREGIONS::iterator lrit = listx.begin(); lrit != listx.end(); lrit++) {
#ifdef _OPENMP
#pragma omp critical
#endif
        { cerr << "running " << lrit->first << "\n"; }
        do_it(lrit->second, lrit->first.data(), b, a, mmqual, mmqual2, pfx);
    }
}
int run_it(LISTOREGIONS listx, char const* a, char const* b, char const* pfx, QUAL_t mmqual, QUAL_t mmqual2) {
    do_it_all(listx, a, b, pfx, mmqual, mmqual2);
    return 0;
}
int run_it_in_fork(LISTOREGIONS& listx, char const* a, char const* b, char const* pfx, QUAL_t mmqual, QUAL_t mmqual2) {
    cerr << "run in fork with backtrace\n";
    struct sigaction act;
    memset(&act, '\0', sizeof(act));
    act.sa_sigaction = sig_action_function;
    act.sa_flags = SA_SIGINFO;
    sigaction(SIGUSR1, &act, 0);
    sigaction(SIGUSR2, &act, 0);
    struct sigaction act2;
    memset(&act2, '\0', sizeof(act2));
    act2.sa_sigaction = crit_err_hdlr;
    act2.sa_flags = SA_RESTART | SA_SIGINFO;
    if (sigaction(SIGSEGV, &act2, (struct sigaction*)NULL) != 0) {
        fprintf(stderr, "error setting signal handler for %d (%s)\n", SIGSEGV, strsignal(SIGSEGV));
        exit(1);
    }
    pid_t child_pid;
    int status = 0;
    while (1) {
        child_pid = fork();
        switch (child_pid) {
        case -1:
            std::cerr << "unable to fork!?! Exiting." << endl;
            exit(1);
            break;
        case 0: {
            cerr << "running within child process " << getpid() << "\n";
            do_it_all(listx, a, b, pfx, mmqual, mmqual2);
            static char messageText[] = "all done";
            union sigval signal_value;
            signal_value.sival_ptr = messageText;
            cerr << "sending shut-down message to parent\n";
            sigqueue(getppid(), SIGUSR2, signal_value);
            sleep(1);
            return 0;
        }
        default:
            std::cerr << "generated child fork pid= " << child_pid << "\n";
            break;
        }
        while ((wait(&status)) > 0) {
            if (status > 0)
                cerr << "SEEMS WE EXCITED IN A MOST UNDIGNIFIED MANNER - SIGNAL : " << status << "\n";
            else {
                std::cerr << "child seems to have finished - should exit here?!?\n";
            }
            if (restart_pos >= sizeof(restart_int) / sizeof(int)) {
                cerr << "will not attempt again\n";
                break;
            }
            cerr << "[" << restart_pos << "] will restart in " << restart_int[restart_pos] << "\n";
            sleep(restart_int[restart_pos]);
            ++restart_pos;
        }
    }
    return 0;
}
string get_vcf_header(char const* SN, int argc, char** argv, bool snp, char const* ref) {
    stringstream hs;
    time_t today = time(0);
    tm* now = localtime(&today);
    char date[9];
    date[8] = '\0';
    sprintf(date, "%d%02d%02d", 1900 + now->tm_year, now->tm_mon, now->tm_mday);
    hs << "##fileformat=VCFv4.0\n##fileDate=" << date << "\n"
                                                         "##source=" << *argv << " " << VERSION << "\n"
                                                                                                   "##command=";
    for (int i = 0; i < argc; ++i) {
        hs << argv[i];
        if (i != argc - 1)
            hs << " ";
    }
    hs << "\n##REFERENCE=" << ref << "\n";
    if (!snp)
        hs << "##INFO=<ID=P,Number=1,Type=Float,Description=\"Indel p value\">\n";
    hs << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
          "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n"
          "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n"
          "##FORMAT=<ID=RR,Number=1,Type=Integer,Description=\"Reference Read Depth\">\n"
          "##FORMAT=<ID=VR,Number=1,Type=Integer,Description=\"Major Variant Read Depth\">\n";
    if (snp)
        hs << "##FILTER=<ID=low_snpqual,Description=\"SNP posterior probability lower than " << ill::SLX_SNP_CUTOFF
           << "\">\n"
              "##FILTER=<ID=equal_majority,Description=\"The called SNP had an equal number of reads indicating "
              "another variant call and base was choosen by highest sumed quality score - this is a deliberate "
              "backwards incompatibility\">\n"
              "##FILTER=<ID=PASS,Description=\"I'll give you 3 guesses...\">\n";
    else
        hs << "##FILTER=<ID=low_qual,Description=\"indel posterior probablitity is less than " << opts::INDEL_PR_CUTOFF
           << "\">\n";
    hs << "##FILTER=<ID=low_VariantReads,Description=\"Number of variant reads is less than "
       << (snp ? opts::SNP_MIN_COV : opts::MIN_VAR_READS)
       << "\">\n"
          "##FILTER=<ID=low_VariantRatio,Description=\"Variant read ratio is less than "
       << (snp ? opts::SNP_HET_MIN : opts::MIN_VAR_RATIO)
       << "\">\n"
          "##FILTER=<ID=low_coverage,Description=\"Total coverage is less than "
       << (snp ? ill::SLX_SNP_MIN_COVERAGE : opts::MIN_DEPTH_COVERAGE) << "\">\n";
    if (snp)
        hs << "##FILTER=<ID=high_coverage,Description=\"Total coverage is more than "
           << opts::SLX_SNP_HIGH_MAP_QUAL_MAX_COVERAGE
           << "\">\n"
              "##FILTER=<ID=single_strand,Description=\"All variant reads are in a single strand direction\">\n";
    else
        hs << "##FILTER=<ID=read_end_ratio,Description=\"Ratio of variants reads within " << (NEAR_END + 2)
           << "bp is greater than " << opts::MAX_NEAR_READ_END_RATIO << "\">\n";
    hs << "##FILTER=<ID=no_data,Description=\"No valid reads on this site\">\n"
          "##FILTER=<ID=no_var,Description=\"No valid variants reads on this site\">\n"
          "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << SN << "\n";
    return hs.str();
}
#define LOGIT(a) cerr << "using " << #a << " of " << a << "\n";
int main(int argc, char** argv) {
    char* a = 0, *b = 0, *c = 0, *pfx = 0, *sample_name = 0;
    QUAL_t mmqual = 0, mmqual2 = 1;
    int tmpmmqual = 0, tmpmmqual2 = 1, threads = 1;
    if (argc == 1) {
        cerr << MY_NAME << ", " << VERSION << ", Daniel S. T. Hughes\n" __FILE__ ", " __DATE__ ", " __TIME__ "\n"
             << "g++-" << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__ << " (" << __cplusplus
             << ")\n\n";
        cerr << "Complete reimplementation, merging, fixes, optimisation, extension, model\nupdate "
                "of AtlasSNP2 & AtlasINDEL2 procedures. See Challis et al., BMC\nBioinformatics, 2012 13:8 & Shen et "
                "al., Genome Res. 2010. 20: 273-280. ";
        cerr << "This\nPre-Release of " << MY_NAME
             << " is available from https://github.com/dsthughes/xatlas.git\n"
                "and released under CREATIVE COMMONS, Attribution-NonCommercial-NoDerivs 3.0\nUnited States. See "
                "https://creativecommons.org/licenses/by-nc-nd/3.0/.\n\n";
        cerr << OPTIONS;
        cerr << "\n";
        exit(0);
    } else if (argc < 4 || argc % 2 == 0)
        std::cerr << "i need args?!?\n", exit(1);
    char options[] = "--";
    bool cleanup = true;
    cerr << "hello, i'm " << *argv << "\n";
    LOGIT(ill::SLX_INDEL_LOW_MAP_QUAL_MAX_COVERAGE);
    LOGIT(ill::SLX_SNP_HIGH_MAP_QUAL_MAX_COVERAGE);
    LOGIT(ill::SLX_SNP_MIN_COVERAGE);
    for (int i = 1; i < argc; i += 2) {
        if (strlen(*(argv + i)) <= sizeof(options) || memcmp(*(argv + i), options, sizeof(options) != 0))
            std::cerr << "this isn't an option " << *(argv + i) << "\n", exit(1);
        if (strcmp(*(argv + i) + sizeof(options) - 1, "ref") == 0)
            a = *(argv + i + 1);
        else if (strcmp(*(argv + i) + sizeof(options) - 1, "bam") == 0)
            b = *(argv + i + 1);
        else if (strcmp(*(argv + i) + sizeof(options) - 1, "scf") == 0)
            c = *(argv + i + 1);
        else if (strcmp(*(argv + i) + sizeof(options) - 1, "samplename") == 0)
            sample_name = *(argv + i + 1);
        else if (strcmp(*(argv + i) + sizeof(options) - 1, "prefix") == 0)
            pfx = *(argv + i + 1);
        else if (strcmp(*(argv + i) + sizeof(options) - 1, "capturebed") == 0)
            opts::capturebed = *(argv + i + 1);
        else if (strcmp(*(argv + i) + sizeof(options) - 1, "depthdebug") == 0 && (strcmp(*(argv + i + 1), "true") == 0))
            opts::depth = true;
        else if (strcmp(*(argv + i) + sizeof(options) - 1, "fork") == 0 && (strcmp(*(argv + i + 1), "true") == 0))
            opts::run_in_fork = true;
        else if (strcmp(*(argv + i) + sizeof(options) - 1, "precalcp") == 0 && (strcmp(*(argv + i + 1), "false") == 0))
            opts::silly_scavenge = false;
        else if (strcmp(*(argv + i) + sizeof(options) - 1, "ignorelowqual") == 0 &&
                 (strcmp(*(argv + i + 1), "false") == 0))
            opts::ignore_low_qual = false;
        else if (strcmp(*(argv + i) + sizeof(options) - 1, "force") == 0 && (strcmp(*(argv + i + 1), "true") == 0))
            opts::force = true;
        else if (strcmp(*(argv + i) + sizeof(options) - 1, "allnonref") == 0 && (strcmp(*(argv + i + 1), "true") == 0))
            opts::allnonref = true;
        else if (strcmp(*(argv + i) + sizeof(options) - 1, "allindels") == 0 && (strcmp(*(argv + i + 1), "true") == 0))
            opts::allindels = true;
        else if (strcmp(*(argv + i) + sizeof(options) - 1, "checkpoint") == 0 &&
                 (strcmp(*(argv + i + 1), "false") == 0))
            opts::checkpoint = false;
        else if (strcmp(*(argv + i) + sizeof(options) - 1, "cleanup") == 0 && (strcmp(*(argv + i + 1), "false") == 0))
            cleanup = false;
        else if (strcmp(*(argv + i) + sizeof(options) - 1, "rate") == 0 && (strcmp(*(argv + i + 1), "true") == 0))
            opts::rate = true;
        else if (strcmp(*(argv + i) + sizeof(options) - 1, "threads") == 0)
            threads = atoi(*(argv + i + 1));
        else if (strcmp(*(argv + i) + sizeof(options) - 1, "minindelmapq") == 0)
            tmpmmqual = atoi(*(argv + i + 1));
        else if (strcmp(*(argv + i) + sizeof(options) - 1, "minsnpmapq") == 0)
            tmpmmqual2 = atoi(*(argv + i + 1));
        else if (strcmp(*(argv + i) + sizeof(options) - 1, "maxcov") == 0)
            opts::SLX_SNP_HIGH_MAP_QUAL_MAX_COVERAGE = atoi(*(argv + i + 1));
        else if (strcmp(*(argv + i) + sizeof(options) - 1, "blockabslim") == 0)
            opts::block_abs_lim = atoi(*(argv + i + 1));
        else if (strcmp(*(argv + i) + sizeof(options) - 1, "blockrellim") == 0)
            opts::block_rel_lim = atof(*(argv + i + 1));
        else if (strcmp(*(argv + i) + sizeof(options) - 1, "nonreffrac") == 0)
            opts::max_alt_frac = atof(*(argv + i + 1));
        else if (strcmp(*(argv + i) + sizeof(options) - 1, "nonreffraccutoff") == 0)
            opts::VRCUTOFF = atoi(*(argv + i + 1));
        else if (strcmp(*(argv + i) + sizeof(options) - 1, "xindel") == 0 && (strcmp(*(argv + i + 1), "true") == 0))
            opts::indeltest = true;
        else if (strcmp(*(argv + i) + sizeof(options) - 1, "gvcftest") == 0 && (strcmp(*(argv + i + 1), "true") == 0))
            opts::gvcftest = true;
        else if (strcmp(*(argv + i) + sizeof(options) - 1, "novar") == 0 && (strcmp(*(argv + i + 1), "true") == 0))
            opts::gvcftestdebug = true;
        else if (strcmp(*(argv + i) + sizeof(options) - 1, "knownindelstest") == 0)
            opts::known_indels = *(argv + i + 1);
        else
            cout << "what is '" << argv[i] << "'?\n", exit(1);
    }
    LOGIT(opts::max_alt_frac);
    LOGIT(opts::block_abs_lim);
    LOGIT(opts::block_rel_lim);
    snprintf(opts::block_label + strlen(opts::block_label), 25, "%dp%da", int(100 * opts::block_rel_lim),
             opts::block_abs_lim);
    LOGIT(opts::block_label);
    LOGIT(opts::SLX_SNP_HIGH_MAP_QUAL_MAX_COVERAGE);
    cout << "let's go\n";
    if (opts::SLX_SNP_HIGH_MAP_QUAL_MAX_COVERAGE > ((uint16_t)-1) >> 3)
        cerr << "max cov too high?!?\n", exit(1);
    if (opts::capturebed && (!opts::gvcftest && !opts::gvcftestdebug))
        cerr << "unless trying to call specific regions alone with WGS etc., capture without gvcf or cf is pretty "
                "silly\n",
            sleep(3);
    if (!opts::capturebed && opts::gvcftestdebug)
        cerr << "don't be silly!?!?\n", exit(1);
    if (opts::gvcftest && opts::gvcftestdebug)
        cerr << "make up your mind\n", exit(1);
    if (a && b && sample_name)
        cerr << "using bam '" << b << "' with ref '" << a << "' "
             << " sample name " << sample_name << "\n";
    else
        cerr << "must set bam, ref, scf & samplename\n", exit(1);
    if (tmpmmqual < 0 || tmpmmqual > 60)
        cerr << "please be serious", exit(1);
    mmqual = tmpmmqual;
    if (tmpmmqual2 < 0 || tmpmmqual2 > 60)
        cerr << "please be serious", exit(1);
    mmqual2 = tmpmmqual2;
    if (tmpmmqual2 < tmpmmqual)
        cerr << "haven't implemented lower mapping qual for snps\n", exit(1);
    if (!pfx)
        cerr << "give me an output prefix\n", exit(1);
    cerr << "using prefix " << pfx << "\n";
    cerr << "using minimum snp map qual " << (short)mmqual << "\n";
    cerr << "using minimum indel map qual " << (short)mmqual2 << "\n";
    cerr << "running with " << threads << " threads\n";
#ifdef _OPENMP
    omp_set_num_threads(threads);
#endif
    LISTOREGIONS list_x;
    std::vector<std::string> bam_scf_order;
    if (c) {
        list_x.insert(make_pair(c, LIST()));
    } else {
        cout << "using bam header for input sequence list\n";
        {
            samFile* in = 0;
            in = sam_open_format(b, "rb", 0);
            assert(in != 0);
            bam_hdr_t* header = sam_hdr_read(in);
            cerr << "there are " << header->n_targets << " scaffolds in bam\n";
            for (int i = 0; i < header->n_targets; ++i) {
                cerr << "scaffold " << i << " : " << header->target_name[i] << " : " << header->target_len[i] << "\n";
                bam_scf_order.push_back(header->target_name[i]);
            }
            bam_hdr_destroy(header);
            sam_close(in);
        }
        if (opts::capturebed) {
            cerr << "using " << opts::capturebed << " for input sequence list\n";
            if (!utils::isregfile(opts::capturebed))
                cerr << "there is a problem with " << opts::capturebed << "\n", exit(1);
            ifstream bf(opts::capturebed);
            if (!bf)
                cerr << "problem reading " << opts::capturebed << "\n";
            string tmp;
            std::string cscf = "";
            while (getline(bf, tmp)) {
                std::vector<std::string> tmpvs;
                tokenise(tmpvs, tmp, '\t');
                if (tmpvs.size() != 3)
                    cerr << "this doesn't look like bed?!?\n", exit(1);
                uint32_t ax = atoi(tmpvs[1].data()), bx = atoi(tmpvs[2].data());
                list_x[tmpvs[0]].push_back(make_pair(ax, bx));
            }
            bf.close();
            std::vector<std::string> tmpscf;
            for (unsigned int i = 0; i < bam_scf_order.size(); ++i)
                if (list_x.count(bam_scf_order[i]) > 0)
                    tmpscf.push_back(bam_scf_order[i]);
            std::swap(bam_scf_order, tmpscf);
        } else
            for (unsigned int i = 0; i < bam_scf_order.size(); ++i)
                list_x.insert(make_pair(bam_scf_order[i], LIST()));
    }
    if (list_x.size() == 0)
        cerr << "unable to identify scaffold " << c << " in bam header\n", exit(1);
    cerr << "will analyse " << list_x.size() << " scaffolds\n";
    if (!opts::run_in_fork)
        run_it(list_x, a, b, pfx, mmqual, mmqual2);
    else
        run_it_in_fork(list_x, a, b, pfx, mmqual, mmqual2);
    cerr << "put in VCF header\n";
    ofstream fsnpv, findelv, fdepth;
    stringstream csnp, cindel, cdepth;
    fsnpv.open((string(pfx) + "_FINAL_snp.vcf").data(), std::fstream::out | std::fstream::trunc);
    findelv.open((string(pfx) + "_FINAL_indel.vcf").data(), std::fstream::out | std::fstream::trunc);
    if (opts::depth) {
        fdepth.open((std::string(pfx) + "_FINAL.depth").data(), std::fstream::out | std::fstream::trunc);
        if (!fdepth)
            cerr << "couldn't open final depth output file\n", exit(1);
        fdepth << "Chr\tPos\tAS_VR\tAS_AR\tAS_RR\tAS_DP\tAID_VR\tAID_AR\tAID_RR\tAID_DP\n";
        fdepth.close();
        cdepth << "cat ";
    }
    if (!fsnpv)
        cerr << "couldn't open final snp output file\n", exit(1);
    if (!findelv)
        cerr << "couldn't open final indel output file\n", exit(1);
    fsnpv << get_vcf_header(sample_name, argc, argv, true, a);
    findelv << get_vcf_header(sample_name, argc, argv, false, a);
    fsnpv.close();
    findelv.close();
    csnp << "cat ";
    cindel << "cat ";
    for (unsigned int i = 0; i < bam_scf_order.size(); ++i) {
        csnp << pfx << "_seq" << bam_scf_order[i] << "_snp.vcf ";
        cindel << pfx << "_seq" << bam_scf_order[i] << "_indel.vcf ";
    }
    csnp << " >> " << pfx << "_FINAL_snp.vcf";
    cindel << " >> " << pfx << "_FINAL_indel.vcf";
    cerr << "concatenating individual scaffolds\n";
    if (system(csnp.str().data()))
        cerr << "unable to generate final output snp file\n";
    if (system(cindel.str().data()))
        cerr << "unable to generate final output indel file\n";
    cerr << "checking file consistency\n";
    if (!check_file_okay((string(pfx) + "_FINAL_snp.vcf").data()) ||
        !check_file_okay((string(pfx) + "_FINAL_indel.vcf").data()))
        cerr << "there's a problem with output file consistency\n", exit(1);
    cerr << "cleaning up\n";
    std::vector<string> cleanuplist;
    stringstream fn;
    if (cleanup) {
        for (LISTOREGIONS::iterator fit = list_x.begin(); fit != list_x.end(); ++fit) {
            fn << pfx << "_seq" << fit->first << "_snp.vcf ";
            cleanuplist.push_back(fn.str());
            fn.str("");
            fn.clear();
            fn << pfx << "_seq" << fit->first << "_indel.vcf ";
            cleanuplist.push_back(fn.str());
            if (opts::depth) {
                fn.str("");
                fn.clear();
                fn << pfx << "_seq" << fit->first << ".depth ";
                cleanuplist.push_back(fn.str());
            }
            if (opts::checkpoint) {
                fn.str("");
                fn.clear();
                fn << pfx << "_seq" << fit->first;
                fn << ".done";
                cleanuplist.push_back(fn.str());
            }
            fn.str("");
            fn.clear();
        }
        fn.str("");
        fn.clear();
        fn << "rm ";
        for (unsigned int i = 0; i < cleanuplist.size(); ++i)
            fn << cleanuplist[i] << " ";
        if (system(fn.str().data()))
            cerr << "unable to cleanup\n";
    }
    return 0;
}
