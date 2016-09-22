#include "EventHolder.hpp"
#include <cmath>

EventHolder::EventHolder()
{
}

EventHolder::EventHolder(pos_t near_end, double snp_max_sub, double snp_max_gap)
    : _near_end(near_end),
      _snp_max_sub(snp_max_sub),
      _snp_max_gap(snp_max_gap),
      last_call_indel(-1),
      last_call_snp(-1)
{
}

EventHolder::EventHolder(const EventHolder &) = default;

EventHolder::~EventHolder() = default;

void EventHolder::addSnp(pos_t pos, const SnpEvent &sv)
{
    this->snps[pos].push_back(sv);
}

void EventHolder::addIndel(pos_t pos, const IndelEvent &iv)
{
    const auto &iv_it = this->indels.find(pos);
    if (iv_it == this->indels.end())
        this->indels[pos].emplace(std::make_pair(iv.id, iv));
    else
    {
        const auto &ivg_it = iv_it->second.find(iv.id);
        if (ivg_it == iv_it->second.end())
            iv_it->second.emplace(std::make_pair(iv.id, iv));
        else
            ivg_it->second.addReadIndelEvent(iv);
    }
}

void EventHolder::collectIndels(bam1_t *read, const Sequence &sequences, CoverageCounter &coverages)
{
    uint32_t var_count = 0, prev_op = UINT32_MAX, *cit = bam_get_cigar(read);
    uint32_t *end_cigar = cit + read->core.n_cigar - 1;
    pos_t indel_offset = 0, indel_pos = read->core.pos, read_pos = 0, ref_pos = 0;
    pos_t leading_sclip = (bam_cigar_op(*cit) == BAM_CSOFT_CLIP)
                              ? bam_cigar_oplen(*cit)
                              : 0;
    pos_t tailing_sclip = (bam_cigar_op(*end_cigar) == BAM_CSOFT_CLIP)
                              ? bam_cigar_oplen(*end_cigar)
                              : 0;
    pos_t rlen = read->core.l_qseq - leading_sclip - tailing_sclip;
    uint8_t *bseq = bam_get_seq(read), *bqual = bam_get_qual(read) + leading_sclip;
    bool strand = (bam_is_rev(read) == 0);
    pos_t c_len, pos1, pos2, end_pos;
    std::string seq;
    AReadsIndelList this_reads_indels;

    while (cit <= end_cigar)
    {
        uint32_t c_op = bam_cigar_op(*cit);

        // coverage and events
        switch (c_op)
        {
        case BAM_CMATCH:
            c_len = (pos_t)bam_cigar_oplen(*cit);
            coverages.addCoverageRange(indel_pos, COVSIDX_INDEL_RR, c_len);
            coverages.addCoverageRange(indel_pos, COVSIDX_INDEL_RR_INS, c_len);
            coverages.addCoverageRange(indel_pos, COVSIDX_INDEL_DP, c_len);
            pos1 = read_pos + leading_sclip;
            pos2 = read->core.pos + ref_pos;
            end_pos = pos1 + c_len;
            while (pos1 < end_pos)
            {
                if (seq_nt16_str[bam_seqi(bseq, pos1)] != sequences.seq[pos2])
                    ++var_count;
                ++pos1;
                ++pos2;
            }
            indel_pos += c_len;
            indel_offset += c_len;
            ref_pos += c_len;
            read_pos += c_len;
            break;
        case BAM_CREF_SKIP:
        case BAM_CEQUAL:
            c_len = bam_cigar_oplen(*cit);
            coverages.addCoverageRange(indel_pos, COVSIDX_INDEL_RR, c_len);
            coverages.addCoverageRange(indel_pos, COVSIDX_INDEL_RR_INS, c_len);
            coverages.addCoverageRange(indel_pos, COVSIDX_INDEL_DP, c_len);
            indel_pos += c_len;
            indel_offset += c_len;
            ref_pos += c_len;
            read_pos += c_len;
            break;
        case BAM_CDIFF:
            c_len = bam_cigar_oplen(*cit);
            coverages.addCoverageRange(indel_pos, COVSIDX_INDEL_RR, c_len);
            coverages.addCoverageRange(indel_pos, COVSIDX_INDEL_RR_INS, c_len);
            coverages.addCoverageRange(indel_pos, COVSIDX_INDEL_DP, c_len);
            var_count += c_len;
            indel_pos += c_len;
            indel_offset += c_len;
            ref_pos += c_len;
            read_pos += c_len;
            break;
        case BAM_CINS:
            c_len = bam_cigar_oplen(*cit);
            pos1 = indel_offset + leading_sclip;
            end_pos = pos1 + c_len;
            seq.clear();
            while (pos1 < end_pos)
            {
                seq.push_back(seq_nt16_str[bam_seqi(bseq, pos1)]);
                ++pos1;
            }
            this_reads_indels.push_back(IndelEvent(indel_pos,
                                                   c_len,
                                                   indel_offset,
                                                   seq,
                                                   strand,
                                                   read->core.qual));
            if (prev_op == UINT32_MAX && indel_pos >= 1)
                coverages.addCoverage(indel_pos - 1, COVSIDX_INDEL_RR_INS);
            ++var_count;
            indel_offset += c_len;
            read_pos += c_len;
            break;
        case BAM_CDEL:
            c_len = bam_cigar_oplen(*cit);
            seq.clear();
            this_reads_indels.push_back(IndelEvent(indel_pos,
                                                   c_len,
                                                   indel_offset,
                                                   seq,
                                                   strand,
                                                   read->core.qual));
            coverages.addCoverageRange(indel_pos, COVSIDX_INDEL_DP, c_len);
            ++var_count;
            indel_pos += c_len;
            ref_pos += c_len;
            break;
        case BAM_CPAD:
            c_len = bam_cigar_oplen(*cit);
            indel_pos += c_len;
            ref_pos += c_len;
            break;
        case BAM_CSOFT_CLIP:
        case BAM_CHARD_CLIP:
            break;
        default:
            break;
        }

        prev_op = c_op;
        ++cit;
    }

    //double var_rate_gap_and_mismatch = round(((double)var_count / rlen) * 10000.0) / 10000.0;
    double var_rate_gap_and_mismatch = (double)var_count / rlen;

    for (auto &iv_it : this_reads_indels)
    {
        end_pos = iv_it.offset;
        if (!iv_it.seq.empty())
            end_pos += iv_it.var_len;

        pos1 = iv_it.offset - this->_near_end - 2;
        if (pos1 < 0)
            pos1 = 0;

        pos2 = end_pos + this->_near_end;
        if (pos2 >= rlen)
            pos2 = rlen - 1;

        big_qual_t qualsum = 0;
        for (ref_pos = pos1; ref_pos <= pos2; ++ref_pos)
            qualsum += bqual[ref_pos];

        iv_it.near_read_end_count = (iv_it.offset + 1 <= this->_near_end || rlen - end_pos < this->_near_end + 2 ? 1 : 0);
        iv_it.avg_nbq = qualsum / (pos2 - pos1 + 1);
        iv_it.var_rate_gap_and_mismatch = var_rate_gap_and_mismatch;
        this->addIndel(iv_it.var_start, iv_it);
    }
}

inline bool EventHolder::snpNqs(pos_t snp_3p, pos_t dist, uint8_t *qquals, pos_t qlen)
{
    bool nqs = false;

    if (qlen > 10)
    {
        if (dist > qlen - 6)
            nqs = true;
        else if (snp_3p >= 5 && qlen - snp_3p >= 6 && qquals[snp_3p] >= 20)
        {
            nqs = true;
            for (int32_t i = -5; i <= 5; ++i)
            {
                if (i != 0 && qquals[snp_3p + i] < 15)
                {
                    nqs = false;
                    break;
                }
            }
        }
    }

    return nqs;
}

void EventHolder::collectSnps(bam1_t *read, const Sequence &sequences, CoverageCounter &coverages)
{
    uint8_t snp_nt_gap = 0;
    int32_t snp_cql = 0;
    uint32_t *cit = bam_get_cigar(read);
    uint32_t *end_cigar = cit + read->core.n_cigar - 1;
    pos_t len = read->core.l_qseq;

    // snp sub and gap
    while (cit <= end_cigar)
    {
        pos_t c_len;

        switch (bam_cigar_op(*cit))
        {
        case BAM_CINS:
            c_len = bam_cigar_oplen(*cit);
            snp_nt_gap += c_len;
            snp_cql += c_len;
            break;
        case BAM_CDEL:
            snp_nt_gap += bam_cigar_oplen(*cit);
            break;
        case BAM_CSOFT_CLIP:
            c_len = bam_cigar_oplen(*cit);
            len -= c_len;
            snp_cql += c_len;
            break;
        case BAM_CMATCH:
            snp_cql += bam_cigar_oplen(*cit);
            break;
        case BAM_CDIFF:
        case BAM_CEQUAL:
        case BAM_CPAD:
        case BAM_CREF_SKIP:
        case BAM_CHARD_CLIP:
            break;
        default:
            break;
        }

        ++cit;
    }

    if (snp_cql != read->core.l_qseq)
    {
        //std::cerr << "CIGAR length/read length mismatch!" << std::endl;
        return;
    }

    uint8_t snp_nt_sub = bam_aux2i(bam_aux_get(read, "NM")) - snp_nt_gap;

    if ((double)snp_nt_sub / len > this->_snp_max_sub ||
        (double)snp_nt_gap / len > this->_snp_max_gap)
    {
        return;
    }

    // coverage and events
    bool strand = (bam_is_rev(read) == 0);
    pos_t c_len, snp_tplace = read->core.pos, snp_qplace = 0;
    cit = bam_get_cigar(read);
    uint8_t snp_nt_snp = 0, *bseq = bam_get_seq(read), *bqual = bam_get_qual(read);
    AReadsSnpList this_reads_snps;

    while (cit <= end_cigar)
    {
        switch (bam_cigar_op(*cit))
        {
        case BAM_CMATCH:
            c_len = bam_cigar_oplen(*cit);
            if (snp_nt_sub == 0)
            {
                coverages.addCoverageRange(snp_tplace, COVSIDX_SNP_RR, c_len);
                coverages.addCoverageRange(snp_tplace, COVSIDX_SNP_DP, c_len);
                snp_tplace += c_len;
                snp_qplace += c_len;
            }
            else
            {
                pos_t pos1 = snp_tplace, len = snp_tplace + c_len;

                while (snp_tplace < len)
                {
                    char refbase = sequences.seq[snp_tplace];
                    char allele = seq_nt16_str[bam_seqi(bseq, snp_qplace)];

                    if (refbase != allele)
                    {
                        ++snp_nt_snp;

                        if (snp_tplace > pos1)
                        {
                            pos_t pos2 = snp_tplace - pos1;
                            coverages.addCoverageRange(pos1, COVSIDX_SNP_RR, pos2);
                            coverages.addCoverageRange(pos1, COVSIDX_SNP_DP, pos2);
                        }
                        pos1 = snp_tplace + 1;

                        if (refbase != 'N' && allele != 'N')
                        {
                            pos_t snp_3p = strand ? read->core.l_qseq - snp_qplace - 1 : snp_qplace;
                            this_reads_snps.push_back(std::make_pair(snp_tplace,
                                                                     SnpEvent(refbase,
                                                                              allele,
                                                                              bqual[snp_qplace],
                                                                              (double)snp_3p / read->core.l_qseq,
                                                                              this->snpNqs(snp_qplace, snp_3p, bqual, read->core.l_qseq),
                                                                              strand)));
                        }
                    }

                    ++snp_tplace;
                    ++snp_qplace;
                }

                if (snp_tplace > pos1)
                {
                    pos_t pos2 = snp_tplace - pos1;
                    coverages.addCoverageRange(pos1, COVSIDX_SNP_RR, pos2);
                    coverages.addCoverageRange(pos1, COVSIDX_SNP_DP, pos2);
                }
            }
            break;
        case BAM_CINS:
        case BAM_CSOFT_CLIP:
            snp_qplace += bam_cigar_oplen(*cit);
            break;
        case BAM_CDEL:
            c_len = bam_cigar_oplen(*cit);
            coverages.addCoverageRange(snp_tplace, COVSIDX_SNP_DP, c_len);
            snp_tplace += c_len;
            break;
        case BAM_CDIFF:
        case BAM_CEQUAL:
        case BAM_CPAD:
        case BAM_CREF_SKIP:
        case BAM_CHARD_CLIP:
            break;
        default:
            break;
        }

        ++cit;
    }

    if (snp_nt_sub != snp_nt_snp)
        return;

    for (const auto &snp_it : this_reads_snps)
        this->addSnp(snp_it.first, snp_it.second);

    /* Set snp type, not required for current Illumina SNP model because coefficient is 0
    for (auto snp_it = this_reads_snps.begin(); snp_it != this_reads_snps.end(); ++snp_it)
    {
        char refbase = snp_it->second.ref_base;
        char allele = snp_it->second.allele_base;

        AReadsSnpList::iterator snp_it2(snp_it);
        ++snp_it2;
        if (snp_it2 != this_reads_snps.end())
        {
            pos_t pos1 = snp_it2->first;
            if (pos1 - snp_it->first == 1)
            {
                if (snp_it->second.type != SNPTYPE_SWAP)
                    snp_it->second.type = SNPTYPE_MNP;
                if (snp_it2->second.type != SNPTYPE_SWAP)
                    snp_it2->second.type = SNPTYPE_MNP;
            }

            if (pos1 - snp_it->first <= 2 &&
                allele == sequences.seq[pos1] &&
                refbase == snp_it2->second.allele_base)
            {
                snp_it->second.type = SNPTYPE_SWAP;
                snp_it2->second.type = SNPTYPE_SWAP;
            }

            ++snp_it2;
            if (snp_it2 != this_reads_snps.end())
            {
                pos1 = snp_it2->first;
                if (pos1 - snp_it->first <= 2 &&
                    allele == sequences.seq[pos1] &&
                    refbase == snp_it2->second.allele_base)
                {
                    snp_it->second.type = SNPTYPE_SWAP;
                    snp_it2->second.type = SNPTYPE_SWAP;
                }
            }
        }

        this->addSnp(snp_it->first, snp_it->second);
    }
    */
}
