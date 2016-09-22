#ifndef _XATLAS_SNPEVENT_H
#define _XATLAS_SNPEVENT_H

#include "Xatlas.hpp"

enum snp_type
{
    SNPTYPE_SNP = 0,
    SNPTYPE_SWAP,
    SNPTYPE_MNP
};
typedef enum snp_type snp_type_e;

class SnpEvent
{
  private:
    SnpEvent();

  public:
    char ref_base;
    char allele_base;
    double rel_pos;
    bool read_strand;
    qual_t qual;
    snp_type_e type;
    double nqs;

    SnpEvent(char ref, char alt, qual_t qual, double rp, bool nq, bool rs);
    SnpEvent(const SnpEvent &);
    ~SnpEvent();
};

/* AReadsSnp */

typedef std::vector< std::pair< pos_t, SnpEvent > > AReadsSnpList;

#endif /* _XATLAS_SNPEVENT_H */
