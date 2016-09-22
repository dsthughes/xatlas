#include "CoverageCounter.hpp"
#include <cstring>

CoverageCounter::CoverageCounter()
    : _covs(nullptr),
      _quals(nullptr),
      _max_cov(0)
{
}

CoverageCounter::CoverageCounter(uint32_t len, coverage_t max_cov)
    : _covs(new struct covs[len]),
      _quals(new struct quals[len]),
      _max_cov(max_cov)
{
    memset(this->_covs, 0, len * sizeof(struct covs));
    memset(this->_quals, 0, len * sizeof(struct quals));
}

CoverageCounter::CoverageCounter(const CoverageCounter &) = default;

CoverageCounter::~CoverageCounter()
{
    delete[] this->_covs;
    delete[] this->_quals;
}

void CoverageCounter::setP(pos_t pos, quals_idx_e idx, float qual)
{
    this->_quals[pos].qual[idx] = qual;
}

float CoverageCounter::getP(pos_t pos, quals_idx_e idx)
{
    return this->_quals[pos].qual[idx];
}

void CoverageCounter::addCoverage(pos_t pos, covs_idx_e idx)
{
    if (this->_covs[pos].cov[idx] < this->_max_cov)
        ++this->_covs[pos].cov[idx];
}

void CoverageCounter::addCoverage(pos_t pos, covs_idx_e idx, coverage_t add_cov)
{
    if (this->_covs[pos].cov[idx] + add_cov <= this->_max_cov)
        this->_covs[pos].cov[idx] += add_cov;
}

void CoverageCounter::addCoverageRange(pos_t pos, covs_idx_e idx, pos_t len)
{
    pos_t p = pos, end = pos + len;

    while (p < end)
    {
        if (this->_covs[p].cov[idx] < this->_max_cov)
            ++this->_covs[p].cov[idx];
        ++p;
    }
}

void CoverageCounter::addCoverageRange(pos_t pos, covs_idx_e idx, pos_t len, coverage_t add_cov)
{
    pos_t p = pos, end = pos + len;

    while (p < end)
    {
        if (this->_covs[p].cov[idx] + add_cov <= this->_max_cov)
            this->_covs[p].cov[idx] += add_cov;
        ++p;
    }
}

void CoverageCounter::setCoverage(pos_t pos, covs_idx_e idx, coverage_t cov)
{
    if (cov <= this->_max_cov)
        this->_covs[pos].cov[idx] = cov;
}

coverage_t CoverageCounter::getCoverage(pos_t pos, covs_idx_e idx)
{
    return this->_covs[pos].cov[idx];
}

bool CoverageCounter::isToxic(pos_t pos, covs_idx_e idx, pos_t l_qseq)
{
    return (this->_covs[pos].cov[idx] >= this->_max_cov &&
            this->_covs[pos + l_qseq / 2].cov[idx] >= this->_max_cov &&
            this->_covs[pos + l_qseq].cov[idx] >= this->_max_cov);
}
