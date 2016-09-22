#ifndef _XATLAS_INDELEVENT_H
#define _XATLAS_INDELEVENT_H

#include "CoverageCounter.hpp"
#include "Sequence.hpp"
#include "Xatlas.hpp"
#include <cstdint>

class IndelEvent
{
  private:
    bool _has_calling_read_in_pos_strand;
    bool _has_calling_read_in_neg_strand;
    big_qual_t _map_qual;

    IndelEvent();

  public:
    pos_t var_start;
    pos_t var_len;
    pos_t offset;
    std::string seq;
    std::string id;
    uint16_t read_count;
    bool isdel;
    uint16_t near_read_end_count;
    big_qual_t avg_nbq;
    double var_rate_gap_and_mismatch;

    IndelEvent(pos_t start,
               pos_t len,
               pos_t offset,
               std::string &seq_,
               bool strand,
               qual_t mqual,
               uint16_t nearend = 0,
               big_qual_t avgnbq = 0,
               double vrate = 0.0);
    IndelEvent(const IndelEvent &);
    ~IndelEvent();
    void addReadIndelEvent(const IndelEvent &iv);
    bool simpleStrandTest();
    double simpleLocalEntropy(Sequence &sequences);
    double getNearReadEndRatio();
    double getMeanAvnqs();
    double getMeanMapq();
    double getMeanVarRate();
};

typedef std::vector< IndelEvent > AReadsIndelList;

#endif /* _XATLAS_INDELEVENT_H */
