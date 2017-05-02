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
// BandedChainAlignment implementation
//

#pragma once

#include <cassert>

#include <pacbio/align/BandedChainAlignment.h>
#include <pacbio/align/internal/BCAlignBlocks.h>

namespace PacBio {
namespace Align {
namespace Internal {

class BandedChainAlignerImpl
{
public:
    BandedChainAlignerImpl(const BandedChainAlignConfig& config);

public:
    BandedChainAlignment Align(const char* target, const size_t targetLen, const char* query,
                               const size_t queryLen,
                               const std::vector<PacBio::Align::Seed>& seeds);

    void StitchCigars(PacBio::Data::Cigar* global, PacBio::Data::Cigar&& local);

private:
    struct Sequences
    {
        const char* target;
        size_t targetLen;
        const char* query;
        size_t queryLen;
    };

    void AlignGapBlock(const PacBio::Align::Seed& nextSeed);
    void AlignGapBlockImpl(const size_t hLength, const size_t vLength);
    void AlignLastGapBlock(void);

    void AlignSeedBlock(const PacBio::Align::Seed& seed);

    void Initialize(const char* target, const size_t targetLen, const char* query,
                    const size_t queryLen);

    std::vector<PacBio::Align::Seed> MergeSeeds(const std::vector<PacBio::Align::Seed>& seeds);

    BandedChainAlignment Result(void);

private:
    const BandedChainAlignConfig& config_;

    StandardGlobalAlignBlock gapBlock_;
    BandedGlobalAlignBlock seedBlock_;
    PacBio::Data::Cigar globalCigar_;
    int64_t globalScore_;
    size_t gapBlockBeginH_;
    size_t gapBlockBeginV_;
    Sequences sequences_;
};

}  // namespace Internal
}  // namespace Align
}  // namespace PacBio
