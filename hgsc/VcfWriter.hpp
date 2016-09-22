#ifndef _XATLAS_VCFWRITER_H
#define _XATLAS_VCFWRITER_H

#include "CoverageCounter.hpp"
#include "EventHolder.hpp"
#include "Xatlas.hpp"
#include <cstdint>
#include <fstream>
#include <string>

class VcfWriter
{
  private:
    CoverageCounter *coverages;
    Sequence *sequences;
    EventHolder *events;
    std::string region;
    const opts_s *opts;
    const ill_s *ill;

    VcfWriter();
    void addField(std::string &str, const char *to_add);

  public:
    VcfWriter(CoverageCounter &_coverages,
              Sequence &_sequences,
              EventHolder &_events,
              const char *_region,
              const opts_s *_opts,
              const ill_s *_ill);
    VcfWriter(const VcfWriter &);
    ~VcfWriter();
    void outputGvcf(pos_t next_var_pos,
                    pos_t last_call_pos,
                    pos_t region_end_pos,
                    std::ofstream &output,
                    quals_idx_e i_quals_,
                    covs_idx_e i_vr_,
                    covs_idx_e i_rr_,
                    covs_idx_e i_dp_);
    void printSnpBuffer(pos_t next_var_pos, const bed_coord_t &seg, std::ofstream &osnpv);
    void printIndelBuffer(pos_t next_var_pos, const bed_coord_t &seg, std::ofstream &oindelv);
};

#endif /* _XATLAS_VCFWRITER_H */
