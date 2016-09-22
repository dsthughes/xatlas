#include "IndelEvent.hpp"
#include <cmath>
#include <cstring>

IndelEvent::IndelEvent()
    : _has_calling_read_in_pos_strand(false),
      _has_calling_read_in_neg_strand(false),
      _map_qual(0),
      var_start(0),
      var_len(0),
      offset(0),
      seq(""),
      id(""),
      read_count(0),
      isdel(false),
      near_read_end_count(0),
      avg_nbq(0),
      var_rate_gap_and_mismatch(0)
{
}

IndelEvent::IndelEvent(pos_t start,
                       pos_t len,
                       pos_t offset,
                       std::string &seq_,
                       bool strand,
                       qual_t mqual,
                       uint16_t nearend,
                       big_qual_t avgnbq,
                       double vrate)
    : _has_calling_read_in_pos_strand(strand),
      _has_calling_read_in_neg_strand(!strand),
      _map_qual(mqual),
      var_start(start),
      var_len(len),
      offset(offset),
      seq(seq_),
      id(seq_.empty() ? std::to_string(len) : seq_),
      read_count(1),
      isdel(seq_.empty()),
      near_read_end_count(nearend),
      avg_nbq(avgnbq),
      var_rate_gap_and_mismatch(vrate)
{
}

IndelEvent::IndelEvent(const IndelEvent &) = default;

IndelEvent::~IndelEvent() = default;

void IndelEvent::addReadIndelEvent(const IndelEvent &iv)
{
    ++this->read_count;
    if (iv._has_calling_read_in_pos_strand)
        this->_has_calling_read_in_pos_strand = true;
    if (iv._has_calling_read_in_neg_strand)
        this->_has_calling_read_in_neg_strand = true;
    this->near_read_end_count += iv.near_read_end_count;
    this->_map_qual += iv._map_qual;
    this->avg_nbq += iv.avg_nbq;
    this->var_rate_gap_and_mismatch += iv.var_rate_gap_and_mismatch;
}

bool IndelEvent::simpleStrandTest()
{
    return (this->_has_calling_read_in_pos_strand && this->_has_calling_read_in_neg_strand);
}

double IndelEvent::simpleLocalEntropy(Sequence &sequences)
{
    pos_t lower, upper, len;

    if (this->isdel)
    {
        lower = this->var_start <= 10 ? 1 : this->var_start - 9;
        upper = this->var_start + this->var_len + 11;
    }
    else
    {
        lower = this->var_start <= 10 ? 1 : this->var_start - 10;
        upper = this->var_start + 12;
    }
    len = upper - lower;

    if (len < this->var_len)
        return 0.0;

    bool match;
    char *q, *r = sequences.seq + lower - 1;
    double e, s;
    uint16_t tmp[44];
    pos_t j, k, i = 0, n = 0, end = len - this->var_len;

    memset(tmp, 0, 44 * sizeof(uint16_t));

    while (i < end)
    {
        bool found = false;
        char *p = r + i;
        j = 0;

        while (j < n)
        {
            match = true;
            q = r + tmp[2 * j];
            k = 0;

            while (k < this->var_len)
            {
                if (p[k] != q[k])
                {
                    match = false;
                    break;
                }

                ++k;
            }

            if (match)
            {
                found = true;
                ++tmp[2 * j + 1];
                break;
            }

            ++j;
        }

        if (!found)
        {
            tmp[2 * n] = i;
            tmp[2 * n + 1] = 1;
            ++n;
        }

        ++i;
    }

    e = 0.0;
    s = 1.0 / (1.0 + end);

    for (j = 0; j < n; ++j)
    {
        double f = s * (double)tmp[2 * j + 1];
        e -= f * log(f);
    }

    //return (double)round(e * 1000.0) / 1000.0;
    return e;
}

double IndelEvent::getNearReadEndRatio()
{
    return (double)this->near_read_end_count / this->read_count;
}

double IndelEvent::getMeanAvnqs()
{
    //return (double)round(((double)this->avg_nbq / this->read_count) * 100.0) / 100.0;
    return (double)this->avg_nbq / this->read_count;
}

double IndelEvent::getMeanMapq()
{
    //return (double)round(((double)this->_map_qual / this->read_count) * 100.0) / 100.0;
    return (double)this->_map_qual / this->read_count;
}

double IndelEvent::getMeanVarRate()
{
    //return (double)round(((double)this->var_rate_gap_and_mismatch / this->read_count) * 100.0) / 100.0;
    return (double)this->var_rate_gap_and_mismatch / this->read_count;
}
