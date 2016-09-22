#ifndef _XATLAS_BAM_H
#define _XATLAS_BAM_H

#include "faidx.h"
#include "sam.h"

class Bam
{
  private:
    bam_hdr_t *_header;
    hts_idx_t *_idx;

    Bam();

  public:
    samFile *sf;
    hts_itr_t *iter;

    explicit Bam(const char *sf_fn);
    Bam(const Bam &);
    ~Bam();
    void setIter(const char *where);
};

#endif /* _XATLAS_BAM_H */
