#ifndef _XATLAS_EVENTHOLDER_H
#define _XATLAS_EVENTHOLDER_H

#include "IndelEvent.hpp"
#include "SnpEvent.hpp"
#include "Xatlas.hpp"
#include "sam.h"

/* TODO convert these to priority queues? */
typedef std::map< pos_t, std::vector< SnpEvent > > SnpMap;
typedef std::map< pos_t, std::map< std::string, IndelEvent > > IndelMap;

class EventHolder
{
  private:
    pos_t _near_end;
    double _snp_max_sub;
    double _snp_max_gap;

    EventHolder();
    bool snpNqs(pos_t snp_3p, pos_t dist, uint8_t *qquals, pos_t qlen);

  public:
    SnpMap snps;
    IndelMap indels;
    pos_t last_call_indel;
    pos_t last_call_snp;

    EventHolder(pos_t near_end, double snp_max_sub, double snp_max_gap);
    EventHolder(const EventHolder &);
    ~EventHolder();
    void addSnp(pos_t pos, const SnpEvent &sv);
    void addIndel(pos_t pos, const IndelEvent &iv);
    void collectIndels(bam1_t *read, const Sequence &sequences, CoverageCounter &coverages);
    void collectSnps(bam1_t *read, const Sequence &sequences, CoverageCounter &coverages);
};

#endif /* _XATLAS_EVENTHOLDER_H */
