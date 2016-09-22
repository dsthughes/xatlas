#ifndef _XATLAS_BLOCKSTATS_H
#define _XATLAS_BLOCKSTATS_H

#include <cstdint>

class BlockStats
{
  private:
    double _low;       // lower bound
    double _high;      // upper bound
    double _abs_lim;   // absolute limit
    double _rel_lim;   // relative limit coefficient
    double _rel_bound; // relative bound size
    double _sum;       // sum of values
    double _sum_sq;    // sum of squared values

    BlockStats();

  public:
    uint32_t k; // block width in number of bases
    double min; // minimum value
    double max; // maximum value

    BlockStats(double block_abs_lim, double block_rel_lim);
    BlockStats(const BlockStats &);
    ~BlockStats();
    bool breakBlock(double value);
    void addValue(double value);
    void resetBlocker();
    double getStdDev();
    double getMean();
};

#endif /* _XATLAS_BLOCKSTATS_H */
