// Author: Derek Barnett

//
// Support for banded alignment of seed chains
//

#pragma once

#include <string>
#include <vector>

#include <pbcopper/align/Seed.h>
#include <pbcopper/data/Cigar.h>

namespace PacBio {
namespace Align {

/// \brief The BandedChainAlignConfig struct provides various parameters used
///        by the BandedChainAlign algorithm.
///
struct BandedChainAlignConfig
{
public:
    float matchScore_;
    float mismatchPenalty_;
    float gapOpenPenalty_;
    float gapExtendPenalty_;
    size_t bandExtend_;

public:
    static BandedChainAlignConfig Default(void);
};

///
/// \brief The BandedChainAlignment struct contains the result of a
///        call to BandedChainAlign.
///
class BandedChainAlignment
{
public:
    BandedChainAlignConfig config_;
    std::string target_;
    std::string query_;
    std::string alignedTarget_;
    std::string alignedQuery_;
    PacBio::Data::Cigar cigar_;

public:
    BandedChainAlignment(void) = default;
    BandedChainAlignment(const BandedChainAlignConfig& config, std::string target,
                         std::string query, PacBio::Data::Cigar cigar);
    BandedChainAlignment(const BandedChainAlignConfig& config, const char* target,
                         const size_t targetLen, const char* query, const size_t queryLen,
                         const PacBio::Data::Cigar& cigar);

    BandedChainAlignment(const BandedChainAlignment&) = default;
    BandedChainAlignment& operator=(const BandedChainAlignment&) = default;

public:
    float Identity(void) const;
    int64_t Score(void) const;
};

///
/// \brief BandedChainAlign
///
///  Peforms banded alignment over a list of seeds.
///
/// \param target     target (reference) sequence
/// \param targetLen  target length
/// \param query      query sequence
/// \param queryLen   query length
/// \param seeds      pre-computed seeds to guide alignment
/// \param config     algorithm parameters
///
/// \return alignment results (pairwise alignment, score, etc)
///
BandedChainAlignment BandedChainAlign(
    const char* target, const size_t targetLen, const char* query, const size_t queryLen,
    const std::vector<PacBio::Align::Seed>& seeds,
    const BandedChainAlignConfig& config = BandedChainAlignConfig::Default());

///
/// \brief BandedChainAlign
///
///  Peforms banded alignment over a list of seeds.
///
///  This is an overloaded method.
///
/// \param target     target (reference) sequence
/// \param query      query sequence
/// \param seeds      pre-computed seeds to guide alignment
/// \param config     algorithm parameters
///
/// \return alignment results (pairwise alignment, score, etc)
///
inline BandedChainAlignment BandedChainAlign(
    const std::string& target, const std::string& query,
    const std::vector<PacBio::Align::Seed>& seeds,
    const BandedChainAlignConfig& config = BandedChainAlignConfig::Default())
{
    return BandedChainAlign(target.c_str(), target.size(), query.c_str(), query.size(), seeds,
                            config);
}

}  // namespace Align
}  // namespace PacBio
