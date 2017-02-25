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
// SIMD local (Smith-Waterman) alignment score
//

#pragma once

#include <string>
#include <vector>

namespace PacBio {
namespace Align {

class LocalAlignment
{
public:
    LocalAlignment(const int32_t targetBegin, const int32_t targetEnd, const int32_t queryBegin,
                   const int32_t queryEnd, const int32_t mismatches, const uint16_t score,
                   const std::vector<uint32_t>& cigar, const std::string& cigarString);

    LocalAlignment(const int32_t targetBegin, const int32_t targetEnd, const int32_t queryBegin,
                   const int32_t queryEnd, const int32_t mismatches, const uint16_t score,
                   std::vector<uint32_t>&& cigar, std::string&& cigarString);

    LocalAlignment(const LocalAlignment&) = delete;
    LocalAlignment(LocalAlignment&&) = default;
    LocalAlignment& operator=(const LocalAlignment&) = delete;
    LocalAlignment& operator=(LocalAlignment&&) = default;
    ~LocalAlignment(void) = default;

public:
    int32_t TargetBegin(void) const { return targetBegin_; }
    int32_t TargetEnd(void) const { return targetEnd_; }
    int32_t QueryBegin(void) const { return queryBegin_; }
    int32_t QueryEnd(void) const { return queryEnd_; }
    int32_t NumMismatches(void) const { return mismatches_; }
    uint16_t Score(void) const { return score_; }
    std::vector<uint32_t> Cigar(void) const { return cigar_; }
    std::string CigarString(void) const { return cigarString_; }

private:
    int32_t targetBegin_;
    int32_t targetEnd_;
    int32_t queryBegin_;
    int32_t queryEnd_;
    int32_t mismatches_;
    uint16_t score_;
    std::vector<uint32_t> cigar_;
    std::string cigarString_;
};

struct LocalAlignConfig
{
public:
    uint8_t MatchScore;
    uint8_t MismatchPenalty;
    uint8_t GapOpenPenalty;
    uint8_t GapExtendPenalty;

public:
    static LocalAlignConfig Default(void);
};

///
/// \brief LocalAlign
///
/// \param target
/// \param query
/// \param config
///
/// \return
///
LocalAlignment LocalAlign(const std::string& target, const std::string& query,
                          const LocalAlignConfig& config = LocalAlignConfig::Default());

///
/// \brief LocalAlign
///
/// \param target
/// \param queries
/// \param config
///
/// \return
///
std::vector<LocalAlignment> LocalAlign(
    const std::string& target, const std::vector<std::string>& queries,
    const LocalAlignConfig& config = LocalAlignConfig::Default());

}  // namespace Align
}  // namespace PacBio
