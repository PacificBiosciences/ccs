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

// Author: David Alexander, Lance Hepler

#pragma once

#include <string>
#include <vector>
#include <memory>

#include <ConsensusCore/Poa/RangeFinder.hpp>

#include <pacbio/ccs/Interval.h>

// some fwd decls
namespace ConsensusCore
{
    class ScoredMutation;
    class PoaGraph;
}

namespace PacBio {
namespace CCS {

class SdpRangeFinder : public ConsensusCore::detail::SdpRangeFinder
{
protected:
    virtual
    ConsensusCore::detail::SdpAnchorVector
    FindAnchors(const std::string& consensusSequence,
                const std::string& readSequence) const final;
};

//
// Details of how a read aligns to the POA graph
//
struct PoaAlignmentSummary
{
    bool ReverseComplementedRead;
    Interval ExtentOnRead;
    Interval ExtentOnConsensus;
    float AlignmentScore;
    float AlignmentIdentity;

    PoaAlignmentSummary()
        : ReverseComplementedRead{false}
        , ExtentOnRead(0, 0)
        , ExtentOnConsensus(0, 0)
        , AlignmentScore{0}
        , AlignmentIdentity{0}
    {}
};

struct PoaAlignmentOptions
{
    bool ClipBegin;
    bool ClipEnd;
};

//
// Partial order aligner with parsimonious memory usage
//
class SparsePoa
{
public:
    // ReadKey: Integer key representing
    //  >= 0: a read in the POA graph
    //    -1: read could not be inserted into POA graph
    using ReadKey = int;


public:
    SparsePoa();
    ~SparsePoa();

    //
    // Add read, which must already have been oriented to be in the
    // "forward" convention
    //
    ReadKey AddRead(std::string readSequence,
                    const PoaAlignmentOptions& alnOptions=PoaAlignmentOptions(),
                    float minScoreToAdd=0);

    //
    // Find better orientation, (fwd or RC) and add as such
    //
    ReadKey OrientAndAddRead(std::string readSequence,
                             const PoaAlignmentOptions& alnOptions=PoaAlignmentOptions(),
                             float minScoreToAdd=0);

    //
    // Walk the POA and get the optimal consensus path
    //
    std::shared_ptr<const ConsensusCore::PoaConsensus>
    FindConsensus(int minCoverage, std::vector<PoaAlignmentSummary>* summaries=NULL) const;

    //
    // Clean up the POA graph, pruning minority paths, to speed up
    // successive AddRead operations.
    //
    void PruneGraph(float minCoverageFraction);

private:
    void repCheck();

private:
    using Path = std::vector<ConsensusCore::PoaGraph::Vertex>;

    ConsensusCore::PoaGraph* graph_;
    std::vector<Path> readPaths_;
    std::vector<bool> reverseComplemented_;
    SdpRangeFinder* rangeFinder_;
};

} // namespace CCS
} // namespace PacBio
