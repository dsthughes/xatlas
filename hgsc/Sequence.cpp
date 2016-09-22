#include "Sequence.hpp"
#include <iostream>
#include <string>

Sequence::Sequence()
    : seq(nullptr),
      len(0)
{
}

Sequence::Sequence(const char *refseq, const char *region)
{
    faidx_t *fai = fai_load(refseq);
    if (fai == nullptr)
    {
        std::cerr << "Failed to load FASTA reference index" << std::endl;
        return;
    }

    int region_len;
    this->seq = fai_fetch(fai, region, &region_len);
    this->len = region_len;

    for (int i = 0; i < region_len; ++i)
        this->seq[i] = toupper(this->seq[i]);

    fai_destroy(fai);
}

Sequence::Sequence(const Sequence &) = default;

Sequence::~Sequence()
{
    free(this->seq);
}

std::string Sequence::getRefSeqRange(pos_t pos, pos_t len)
{
    return std::string(this->seq + pos, len);
}
