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

// Authors: David Alexander, Lance Hepler

#pragma once

#include <string>
#include <vector>

#include <pacbio/consensus/align/AlignConfig.h>

namespace PacBio {
namespace Consensus {
namespace {

// Utility functions common to implementations of aligners

inline int Max3(int a, int b, int c) { return std::max((a), std::max((b), (c))); }
inline int ArgMax3(int a, int b, int c)
{
    if (a >= b && a >= c)
        return 0;
    else if (b >= c)
        return 1;
    else
        return 2;
}

}  // anonymous namespace

/// \brief A pairwise alignment
class PairwiseAlignment
{
private:
    std::string target_;
    std::string query_;
    std::string transcript_;

public:
    // target string, including gaps; usually the "reference"
    std::string Target() const;

    // query string, including gaps; usually the "read"
    std::string Query() const;

    // transcript as defined by Gusfield pg 215.
    std::string Transcript() const;

public:
    float Accuracy() const;
    int Matches() const;
    int Errors() const;
    int Mismatches() const;
    int Insertions() const;
    int Deletions() const;
    int Length() const;

public:
    PairwiseAlignment(const std::string& target, const std::string& query);

    static PairwiseAlignment* FromTranscript(const std::string& transcript,
                                             const std::string& unalnTarget,
                                             const std::string& unalnQuery);
};

PairwiseAlignment* Align(const std::string& target, const std::string& query, int* score,
                         AlignConfig config = AlignConfig::Default());

PairwiseAlignment* Align(const std::string& target, const std::string& query,
                         AlignConfig config = AlignConfig::Default());

// These calls return an array, same len as target, containing indices into the query string.
std::vector<int> TargetToQueryPositions(const std::string& transcript);
std::vector<int> TargetToQueryPositions(const PairwiseAlignment& aln);

}  // namespace Consensus
}  // namespace PacBio
