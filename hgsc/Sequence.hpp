#ifndef _XATLAS_SEQUENCE_H
#define _XATLAS_SEQUENCE_H

#include "Xatlas.hpp"
#include "faidx.h"

class Sequence
{
  private:
    Sequence();

  public:
    char *seq;
    uint32_t len;

    Sequence(const char *refseq, const char *region);
    Sequence(const Sequence &);
    ~Sequence();
    std::string getRefSeqRange(pos_t pos, pos_t len);
};

#endif /* _XATLAS_SEQUENCE_H */
