// Copyright (c) 2016, Pacific Biosciences of California, Inc.
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
    BandedChainAlignment(const BandedChainAlignConfig& config, const std::string& target,
                         const std::string& query, const PacBio::Data::Cigar& cigar);
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
