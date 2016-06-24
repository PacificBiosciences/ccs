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

// Author: David Alexander

#pragma once

#include <boost/utility.hpp>
#include <climits>
#include <string>
#include <utility>
#include <vector>

#include <pacbio/consensus/align/AlignConfig.h>
#include <pacbio/consensus/poa/PoaGraph.h>

namespace PacBio {
namespace Consensus {

class PoaGraph;
class PoaGraphPath;
class ScoredMutation;

AlignConfig DefaultPoaConfig(AlignMode mode = AlignMode::GLOBAL);

/// \brief A multi-sequence consensus obtained from a partial-order alignment
struct PoaConsensus : private boost::noncopyable
{
    const std::string Sequence;
    PoaGraph Graph;
    std::vector<PoaGraph::Vertex> Path;

    PoaConsensus(const std::string& css, const PoaGraph& g,
                 const std::vector<PoaGraph::Vertex>& ConsensusPath);

    // NB: this constructor exists to provide a means to avoid an unnecessary
    // copy of the
    // boost graph.  If we had move semantics (C++11) we would be able to get by
    // without
    // this.
    PoaConsensus(const std::string& css, const detail::PoaGraphImpl& g,
                 const std::vector<PoaGraph::Vertex>& ConsensusPath);

    ~PoaConsensus();

    static const PoaConsensus* FindConsensus(const std::vector<std::string>& reads);

    static const PoaConsensus* FindConsensus(const std::vector<std::string>& reads,
                                             const AlignConfig& config, int minCoverage = -INT_MAX);

    static const PoaConsensus* FindConsensus(const std::vector<std::string>& reads, AlignMode mode,
                                             int minCoverage = -INT_MAX);

public:
    // Additional accessors, which do things on the graph/graphImpl
    // LikelyVariants

public:
    std::string ToGraphViz(int flags = 0) const;

    void WriteGraphVizFile(std::string filename, int flags = 0) const;
};

}  // namespace Consensus
}  // namespace PacBio
