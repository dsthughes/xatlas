#include "Bam.hpp"
#include <iostream>

Bam::Bam()
    : _header(nullptr),
      _idx(nullptr),
      sf(nullptr),
      iter(nullptr)
{
}

Bam::Bam(const char *sf_fn)
{
    if ((this->sf = sam_open_format(sf_fn, "rb", 0)) == nullptr)
    {
        std::cerr << "Failed to load BAM file \"" << sf_fn << "\"" << std::endl;
        return;
    }

    if ((this->_header = sam_hdr_read(this->sf)) == nullptr)
    {
        std::cerr << "Failed to load header for BAM file \"" << sf_fn << "\"" << std::endl;
        return;
    }

    if ((this->_idx = sam_index_load(this->sf, sf_fn)) == nullptr)
    {
        std::cerr << "Failed to load index for BAM file \"" << sf_fn << "\"" << std::endl;
        return;
    }

    this->iter = nullptr;
}

Bam::Bam(const Bam &) = default;

Bam::~Bam()
{
    hts_itr_destroy(this->iter);
    hts_idx_destroy(this->_idx);
    bam_hdr_destroy(this->_header);
    sam_close(this->sf);
}

void Bam::setIter(const char *where)
{
    this->iter = sam_itr_querys(this->_idx, this->_header, where);
}
