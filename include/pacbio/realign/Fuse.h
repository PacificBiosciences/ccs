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

// Author: Armin Töpfer

#pragma once

#include <fstream>
#include <string>
#include <vector>

#include <pacbio/data/MSA.h>
#include <pbbam/BamRecord.h>

namespace PacBio {
namespace Realign {

class Fuse
{
public:
    Fuse(const std::string& ccsInput);
    Fuse(const std::vector<Data::ArrayRead>& arrayReads);

public:
    std::string ConsensusSequence() const { return consensusSequence_; }

private:
    std::vector<Data::ArrayRead> FetchAlignedReads(const std::string& ccsInput) const;
    std::string CreateConsensus(const std::vector<Data::ArrayRead>& arrayReads) const;
    std::map<int, std::pair<std::string, int>> CollectInsertions(const Data::MSA& msa) const;
    std::pair<int, std::string> FindInsertions(
        std::map<int, std::pair<std::string, int>>* posInsCov, int windowSize = 20) const;

private:
    static constexpr int minCoverage_ = 50;
    static constexpr int minInsertionCoverage_ = 50;

    std::string consensusSequence_;
};
}
}  // ::PacBio::Realign