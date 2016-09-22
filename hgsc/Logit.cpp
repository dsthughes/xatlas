#include "Logit.hpp"
#include <cmath>

/* clang-format off */
double snp_logit(const ill_s *ill,
                 double q,
                 double s,
                 double n,
                 double r)
{
    return 1.0 - (1.0 / (1.0 + exp(ill->SLX_SNP_intercept +
                                   ill->SLX_SNP_quality_score * q +
                                   ill->SLX_SNP_NQS_pass      * s +
                                   ill->SLX_SNP_swap          * n +
                                   ill->SLX_SNP_rel_pos       * r)));
}

double indel_logit_wgs(const ill_s *ill,
                       double local_ref_simple_entropy,
                       double strand_dir,
                       double mean_avnqs,
                       double mean_mapq,
                       double mean_var_rate,
                       coverage_t read_count,
                       coverage_t total_depth)
{
    return 1.0 - (1.0 / (1.0 + exp(ill->WGS_INDEL_intercept +
                                   ill->WGS_INDEL_simple_local_entropy * local_ref_simple_entropy +
                                   ill->WGS_INDEL_strand_dir           * strand_dir +
                                   ill->WGS_INDEL_mean_avg_nqs         * mean_avnqs +
                                   ill->WGS_INDEL_norm_var_square      * ((double)(read_count * read_count) / total_depth) +
                                   ill->WGS_INDEL_mean_map_qual        * mean_mapq +
                                   ill->WGS_INDEL_mean_var_rate        * mean_var_rate)));
}

double indel_logit(const ill_s *ill,
                   double local_ref_simple_entropy,
                   double strand_dir,
                   double mean_avnqs,
                   coverage_t read_count,
                   coverage_t total_depth)
{
    return 1.0 - (1.0 / (1.0 + exp(ill->SLX_INDEL_intercept +
                                   ill->SLX_INDEL_local_entropy * local_ref_simple_entropy +
                                   ill->SLX_INDEL_strand_dir    * strand_dir +
                                   ill->SLX_INDEL_mean_avg_nqs  * mean_avnqs +
                                   ill->SLX_INDEL_norm_var_sq   * ((double)(read_count * read_count) / total_depth))));
}
/* clang-format on */
