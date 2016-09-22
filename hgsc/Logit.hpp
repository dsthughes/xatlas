#ifndef _XATLAS_LOGIT_H
#define _XATLAS_LOGIT_H

#include "Xatlas.hpp"

double snp_logit(const ill_s *ill, double q, double s, double n, double r);
double indel_logit_wgs(const ill_s *ill, double local_ref_simple_entropy, double strand_dir, double mean_avnqs, double mean_mapq, double mean_var_rate, coverage_t read_count, coverage_t total_depth);
double indel_logit(const ill_s *ill, double local_ref_simple_entropy, double strand_dir, double mean_avnqs, coverage_t read_count, coverage_t total_depth);

#endif /* _XATLAS_LOGIT_H */
