// Author: Derek Barnett

//
// SIMD local (Smith-Waterman) alignment
//

#include <iostream>
#include <utility>

#include <ssw_cpp.h>

#include <pacbio/align/LocalAlignment.h>

namespace PacBio {
namespace Align {

static inline LocalAlignment FromSSW(StripedSmithWaterman::Alignment&& sswAl)
{
    return LocalAlignment{
        sswAl.ref_begin,  sswAl.ref_end,  sswAl.query_begin,      sswAl.query_end,
        sswAl.mismatches, sswAl.sw_score, std::move(sswAl.cigar), std::move(sswAl.cigar_string)};
}

LocalAlignConfig LocalAlignConfig::Default() { return LocalAlignConfig{2, 2, 3, 1}; }

LocalAlignment::LocalAlignment(const int32_t targetBegin, const int32_t targetEnd,
                               const int32_t queryBegin, const int32_t queryEnd,
                               const int32_t mismatches, const uint16_t score,
                               std::vector<uint32_t> cigar, std::string cigarString)
    : targetBegin_(targetBegin)
    , targetEnd_(targetEnd)
    , queryBegin_(queryBegin)
    , queryEnd_(queryEnd)
    , mismatches_(mismatches)
    , score_(score)
    , cigar_(std::move(cigar))
    , cigarString_(std::move(cigarString))
{
}

LocalAlignment LocalAlign(const std::string& target, const std::string& query,
                          const LocalAlignConfig& config)
{
    StripedSmithWaterman::Aligner aligner{config.MatchScore, config.MismatchPenalty,
                                          config.GapOpenPenalty, config.GapExtendPenalty};
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment;

    aligner.Align(query.c_str(), target.c_str(), target.size(), filter, &alignment);
    return FromSSW(std::move(alignment));
}

std::vector<LocalAlignment> LocalAlign(const std::string& target,
                                       const std::vector<std::string>& queries,
                                       const LocalAlignConfig& config)
{
    StripedSmithWaterman::Aligner aligner{config.MatchScore, config.MismatchPenalty,
                                          config.GapOpenPenalty, config.GapExtendPenalty};
    StripedSmithWaterman::Filter filter;
    aligner.SetReferenceSequence(target.c_str(), target.size());

    std::vector<LocalAlignment> results;
    results.reserve(queries.size());
    for (const auto& query : queries) {
        StripedSmithWaterman::Alignment alignment;
        aligner.Align(query.c_str(), filter, &alignment);
        results.push_back(FromSSW(std::move(alignment)));
    }
    return results;
}

}  // namespace Align
}  // namespace PacBio
