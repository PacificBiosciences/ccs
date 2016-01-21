// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
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

// Author: Brett Bowman

#pragma once

namespace PacBio {
namespace SparseAlignment {

// Default values for struct fields
static const size_t NUM_CANDIDATES    = 10;
static const long   MIN_SCORE         = 18;
static const int    MATCH_SCORE       =  5;
static const int    NON_MATCH_PENALTY =  0;
static const int    INSERTION_PENALTY = -4;
static const int    DELETION_PENALTY  = -8;
static const int    MAX_SEED_GAP      = 200;

/// A simple struct for represting a complete specialization for an
///  implementation of the Baker-Giancarlo SDP algorithm.  Templating
///  around this struct allows us to greatly reduce the number of 
///  parameters that need to be thrown around.
///
struct ChainSeedsConfig
{
    // Default constructor
    ChainSeedsConfig(const size_t numCandidatesArg,
                     const long minScoreArg,
                     const int matchScoreArg,
                     const int nonMatchPenaltyArg,
                     const int insertionPenaltyArg,
                     const int deletionPenaltyArg,
                     const int maxSeedGapArg)
        : numCandidates(numCandidatesArg)
        , minScore(minScoreArg)
        , matchScore(matchScoreArg)
        , nonMatchPenalty(nonMatchPenaltyArg)
        , insertionPenalty(insertionPenaltyArg)
        , deletionPenalty(deletionPenaltyArg)
        , maxSeedGap(maxSeedGapArg)
        {}

    // Abbreviated constructors
    ChainSeedsConfig(const size_t numCandidatesArg)
        : numCandidates(numCandidatesArg)
        , minScore(MIN_SCORE)
        , matchScore(MATCH_SCORE)
        , nonMatchPenalty(NON_MATCH_PENALTY)
        , insertionPenalty(INSERTION_PENALTY)
        , deletionPenalty(DELETION_PENALTY)
        , maxSeedGap(MAX_SEED_GAP)
        {}

    ChainSeedsConfig()
        : numCandidates(NUM_CANDIDATES)
        , minScore(MIN_SCORE)
        , matchScore(MATCH_SCORE)
        , nonMatchPenalty(NON_MATCH_PENALTY)
        , insertionPenalty(INSERTION_PENALTY)
        , deletionPenalty(DELETION_PENALTY)
        , maxSeedGap(MAX_SEED_GAP)
        {}

    size_t numCandidates;
    long minScore;
    int matchScore;
    int nonMatchPenalty;
    int insertionPenalty;
    int deletionPenalty;
    int maxSeedGap;
};

}}  // PacBio::SparseAlignment
