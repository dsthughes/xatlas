#ifndef _XATLAS_COVERAGECOUNTER_H
#define _XATLAS_COVERAGECOUNTER_H

#include "Xatlas.hpp"

enum covs_idx
{
    COVSIDX_INDEL_VR = 0,
    COVSIDX_INDEL_RR,
    COVSIDX_INDEL_RR_INS,
    COVSIDX_INDEL_AR,
    COVSIDX_INDEL_DP,
    COVSIDX_SNP_VR,
    COVSIDX_SNP_RR,
    COVSIDX_SNP_AR,
    COVSIDX_SNP_DP
};
typedef covs_idx covs_idx_e;

enum quals_idx
{
    QUALSIDX_SNP_P = 0,
    QUALSIDX_INDEL_P
};
typedef quals_idx quals_idx_e;

class CoverageCounter
{
  private:
    struct covs
    {
        uint16_t cov[9];
    };
    struct covs *_covs;

    struct quals
    {
        float qual[2];
    };
    struct quals *_quals;

    coverage_t _max_cov;

    CoverageCounter();

  public:
    CoverageCounter(uint32_t len, coverage_t max_cov);
    CoverageCounter(const CoverageCounter &);
    ~CoverageCounter();
    void setP(pos_t pos, quals_idx_e idx, float qual);
    float getP(pos_t pos, quals_idx_e idx);
    void addCoverage(pos_t pos, covs_idx_e idx);
    void addCoverage(pos_t pos, covs_idx_e idx, coverage_t add_cov);
    void addCoverageRange(pos_t pos, covs_idx_e idx, pos_t len);
    void addCoverageRange(pos_t pos, covs_idx_e idx, pos_t len, coverage_t add_cov);
    void setCoverage(pos_t pos, covs_idx_e idx, coverage_t cov);
    coverage_t getCoverage(pos_t pos, covs_idx_e idx);
    bool isToxic(pos_t pos, covs_idx_e idx, pos_t l_qseq);
};

#endif /* _XATLAS_COVERAGECOUNTER_H */
