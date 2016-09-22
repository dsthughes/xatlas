#include "SnpEvent.hpp"

SnpEvent::SnpEvent()
    : ref_base(0),
      allele_base(0),
      rel_pos(0.0),
      read_strand(false),
      qual(0),
      type(SNPTYPE_SNP),
      nqs(0.0)
{
}

SnpEvent::SnpEvent(char ref,
                   char alt,
                   qual_t qual,
                   double rp,
                   bool nq,
                   bool rs)
    : ref_base(ref),
      allele_base(alt),
      rel_pos(rp),
      read_strand(rs),
      qual(qual),
      type(SNPTYPE_SNP),
      nqs(nq ? 1.0 : 0.0)
{
}

SnpEvent::SnpEvent(const SnpEvent &) = default;

SnpEvent::~SnpEvent() = default;
