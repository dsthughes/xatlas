#include "BlockStats.hpp"
#include <algorithm>
#include <cmath>

BlockStats::BlockStats()
    : _abs_lim(0.0),
      _rel_lim(0.0)
{
    this->resetBlocker();
}

BlockStats::BlockStats(double block_abs_lim, double block_rel_lim)
    : _abs_lim(block_abs_lim),
      _rel_lim(block_rel_lim)
{
    this->resetBlocker();
}

BlockStats::BlockStats(const BlockStats &) = default;

BlockStats::~BlockStats() = default;

double BlockStats::getMean()
{
    return (this->k > 0)
               ? this->_sum / this->k
               : 0.0;
}

/* sample std dev */
double BlockStats::getStdDev()
{
    return (this->k > 1)
               ? sqrt((this->_sum_sq - this->_sum * this->_sum / this->k) / (this->k - 1.0))
               : 0.0;
}

/* would value invalidate current block if added */
bool BlockStats::breakBlock(double value)
{
    return (this->k > 0 &&
            (value < this->_low ||
             value > this->_high ||
             (value < this->min && value + std::max(value * this->_rel_lim, this->_abs_lim) < this->max)));
}

/* value must be valid to extend block */
void BlockStats::addValue(double value)
{
    if (this->k == 0)
    {
        this->min = this->max = value;
        this->_rel_bound = this->_rel_lim * value;

        if (this->_rel_bound > this->_abs_lim)
        {
            this->_low = value - this->_rel_bound;
            this->_high = value + this->_rel_bound;
        }
        else
        {
            this->_low = value - this->_abs_lim;
            this->_high = value + this->_abs_lim;
        }
    }
    else if (value < this->min)
    {
        this->min = value;
        this->_rel_bound = this->_rel_lim * value;

        if (this->_rel_bound > this->_abs_lim)
        {
            this->_low = this->max - this->_rel_bound;
            this->_high = value + this->_rel_bound;
        }
        else
        {
            this->_low = this->max - this->_abs_lim;
            this->_high = value + this->_abs_lim;
        }
    }
    else if (value > this->max)
    {
        this->max = value;
        this->_low = value + this->min - this->_high;
    }

    this->_sum += value;
    this->_sum_sq += value * value;
    ++this->k;
}

void BlockStats::resetBlocker()
{
    this->k = 0;
    this->min = 0.0;
    this->max = 0.0;
    this->_low = 0.0;
    this->_high = 0.0;
    this->_rel_bound = 0.0;
    this->_sum = 0.0;
    this->_sum_sq = 0.0;
}
