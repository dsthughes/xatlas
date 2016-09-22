#ifndef _XATLAS_H
#define _XATLAS_H

#include <cstdint>
#include <map>
#include <string>
#include <vector>

typedef uint16_t coverage_t;
typedef uint8_t qual_t;
typedef uint16_t big_qual_t;
typedef int32_t ref_id_t;
typedef int32_t pos_t;
typedef std::pair< pos_t, pos_t > bed_coord_t;
typedef std::vector< bed_coord_t > bed_coord_list_t;
typedef std::map< std::string, bed_coord_list_t > bed_coord_map_t;
typedef std::pair< std::string, bed_coord_list_t > bed_chr_t;
typedef std::vector< bed_chr_t > bed_list_t;
typedef std::vector< std::pair< std::string, uint32_t > > contigs_list_t;

typedef struct opts
{
    bool indeltest;
    bool allnonref;
    bool gvcftest;
    bool scavenge;
    int32_t NEAR_END;
    coverage_t MIN_VAR_READS;
    bool STRAND_DIR_FILTER;
    coverage_t MIN_DEPTH_COVERAGE;
    double MIN_VAR_RATIO;
    double SNP_MIN_PR;
    double INDEL_PR_CUTOFF;
    double MAX_NEAR_READ_END_RATIO;
    coverage_t INDEL_DEPTH_CUTOFF;
    double INDEL_HOM_VAR_CUTOFF;
    double INDEL_HET_CUTOFF;
    double SNP_HET_MIN;
    double SNP_HET_MAX;
    double SNP_STRAND_RATIO_CUTOFF;
    coverage_t SNP_STRAND_TEST_COV_CUTOFF;
    coverage_t SNP_MIN_COV;
    unsigned block_abs_lim;
    double block_rel_lim;
    double block_rel_min;
    double max_alt_frac;
    pos_t last_call_max;
    coverage_t VRCUTOFF;
    char block_label[24];
    const char *capturebed;
    coverage_t SLX_SNP_HIGH_MAP_QUAL_MAX_COVERAGE;
} opts_s;

typedef struct ill
{
    double SLX_SNP_intercept;
    double SLX_SNP_quality_score;
    double SLX_SNP_NQS_pass;
    double SLX_SNP_swap;
    double SLX_SNP_rel_pos;
    double SLX_SNP_MAX_SUB;
    double SLX_SNP_MAX_INDEL;
    coverage_t SLX_INDEL_LOW_MAP_QUAL_MAX_COVERAGE;
    //coverage_t SLX_SNP_HIGH_MAP_QUAL_MAX_COVERAGE;
    coverage_t SLX_SNP_MIN_COVERAGE;
    double prior_err_c;
    double prior_snp_c;
    double SLX_SNP_CUTOFF;
    double SLX_err_prior_arr[10];
    double SLX_snp_prior_arr[10];

    double SLX_INDEL_intercept;
    double SLX_INDEL_local_entropy;
    double SLX_INDEL_strand_dir;
    double SLX_INDEL_norm_var_sq;
    double SLX_INDEL_mean_avg_nqs;

    double WGS_INDEL_intercept;
    double WGS_INDEL_simple_local_entropy;
    double WGS_INDEL_strand_dir;
    double WGS_INDEL_norm_var_square;
    double WGS_INDEL_mean_avg_nqs;
    double WGS_INDEL_mean_map_qual;
    double WGS_INDEL_mean_var_rate;
} ill_s;

#define LOGIT(a) std::cerr << "Using " << #a << " of " << (a) << std::endl;

#endif /* _XATLAS_H */
