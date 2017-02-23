// Copyright (c) 2011-2016, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: Derek Barnett

//
// SIMD local (Smith-Waterman) alignment
//

#include <iostream>

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

LocalAlignConfig LocalAlignConfig::Default(void) { return LocalAlignConfig{2, 2, 3, 1}; }

LocalAlignment::LocalAlignment(const int32_t targetBegin, const int32_t targetEnd,
                               const int32_t queryBegin, const int32_t queryEnd,
                               const int32_t mismatches, const uint16_t score,
                               const std::vector<uint32_t>& cigar, const std::string& cigarString)
    : targetBegin_(targetBegin)
    , targetEnd_(targetEnd)
    , queryBegin_(queryBegin)
    , queryEnd_(queryEnd)
    , mismatches_(mismatches)
    , score_(score)
    , cigar_(cigar)
    , cigarString_(cigarString)
{
}

LocalAlignment::LocalAlignment(const int32_t targetBegin, const int32_t targetEnd,
                               const int32_t queryBegin, const int32_t queryEnd,
                               const int32_t mismatches, const uint16_t score,
                               std::vector<uint32_t>&& cigar, std::string&& cigarString)
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
