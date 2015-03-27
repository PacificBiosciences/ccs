// Copyright (c) 2014, Pacific Biosciences of California, Inc.
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

#include <stdbool.h>
#include <string>

namespace PBSeqAnalysis {
namespace PBChimera {

/**
 * @brief Label that annotates a read for a single chimeric breakpoint
 */
struct ChimeraLabel {
    // Instance variables
    std::string sequenceId;
    bool chimeraFlag;
    std::string leftParentId;
    std::string rightParentId;
    int32_t crossover;
    double score;

    // Default Constructor
    ChimeraLabel(std::string sequenceIdArg,
                 std::string leftParentArg,
                 std::string rightParentArg,
                 int32_t crossoverArg,
                 double scoreArg)
            : sequenceId(sequenceIdArg)
            , chimeraFlag(false)
            , leftParentId(leftParentArg)
            , rightParentId(rightParentArg)
            , crossover(crossoverArg)
            , score(scoreArg)
            {};

    // Name-Only or Place-Holder Constructor
    explicit ChimeraLabel(std::string sequenceIdArg)
            : sequenceId(sequenceIdArg)
            , chimeraFlag(false)
            , leftParentId("N/A")
            , rightParentId("N/A")
            , crossover(-1)
            , score(-1.0f)
            {};

    // Empty or Dummy Constructor
    ChimeraLabel()
            : sequenceId("Dummy")
            , chimeraFlag(false)
            , leftParentId("N/A")
            , rightParentId("N/A")
            , crossover(-1)
            , score(-1.0f)
            {};

    // Move constructor
    ChimeraLabel(ChimeraLabel&& src) = default;
    // Copy constructor is deleted!
    ChimeraLabel(const ChimeraLabel& src) = default;
    // Move assignment constructor
    ChimeraLabel& operator=(ChimeraLabel&& rhs) = default;
    // Copy assignment constructor is deleted!
    ChimeraLabel& operator=(const ChimeraLabel& rhs) = default;
    // Destructor
    ~ChimeraLabel() = default;
};

std::ostream& operator<<(std::ostream &o, const ChimeraLabel& label)
{
    // Stream the Sequence Id first
    o << label.sequenceId << ",";

    // Then a human-readable representation of the flag
    if (label.chimeraFlag)
        o << "True"       << ",";
    else
        o << "False"      << ",";

    // The score is only meaningfully defined > 0
    if (label.score > 0.0f)
        o << label.score  << ",";
    else
        o << "NaN"        << ",";

    // Finally the parents and the putative crossover
    o << label.leftParentId  << ","
      << label.rightParentId << ","
      << label.crossover;

    // Return the stream reference
    return o;
}

}  // namespace PBChimera
}  // namespace PBSeqAnalysis
