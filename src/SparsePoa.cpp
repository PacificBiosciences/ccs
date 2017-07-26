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

// Author: Lance Hepler

#include <algorithm>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <pacbio/align/AlignConfig.h>
#include <pacbio/data/Sequence.h>
#include <pacbio/denovo/PoaConsensus.h>
#include <pacbio/denovo/PoaGraph.h>

#include <pacbio/ccs/SparseAlignment.h>
#include <pacbio/denovo/SparsePoa.h>
#include <pbcopper/logging/Logging.h>

using PacBio::Poa::detail::SdpAnchorVector;
using PacBio::Align::AlignConfig;
using PacBio::Align::AlignMode;
using PacBio::Poa::DefaultPoaConfig;
using PacBio::Poa::PoaConsensus;
using PacBio::Poa::PoaGraph;
using PacBio::Data::ReverseComplement;

namespace PacBio {
namespace Poa {

typedef PoaGraph::Vertex Vertex;

SdpAnchorVector SdpRangeFinder::FindAnchors(const std::string& consensusSequence,
                                            const std::string& readSequence) const
{
    return CCS::SparseAlign(6, consensusSequence, readSequence);
}

SparsePoa::SparsePoa()
    : graph_(new PoaGraph())
    , readPaths_()
    , reverseComplemented_()
    , rangeFinder_(new SdpRangeFinder())
{
}

SparsePoa::~SparsePoa()
{
    delete graph_;
    delete rangeFinder_;
}

SparsePoa::ReadKey SparsePoa::AddRead(const std::string& readSequence,
                                      const PoaAlignmentOptions& /* alnOptions */,
                                      float minScoreToAdd)
{
    AlignConfig config = DefaultPoaConfig(AlignMode::LOCAL);
    Path outputPath;
    ReadKey key = -1;

    if (graph_->NumReads() == 0) {
        graph_->AddFirstRead(readSequence, &outputPath);
        readPaths_.push_back(outputPath);
        reverseComplemented_.push_back(false);
        key = graph_->NumReads() - 1;
    } else {
        auto c = graph_->TryAddRead(readSequence, config, rangeFinder_);

        if (c->Score() >= minScoreToAdd) {
            graph_->CommitAdd(c, &outputPath);
            readPaths_.push_back(outputPath);
            reverseComplemented_.push_back(false);
            key = graph_->NumReads() - 1;
        }

        delete c;
    }

    return key;
}

SparsePoa::ReadKey SparsePoa::OrientAndAddRead(const std::string& readSequence,
                                               const PoaAlignmentOptions& /* alnOptions */,
                                               float minScoreToAdd)
{
    AlignConfig config = DefaultPoaConfig(AlignMode::LOCAL);
    Path outputPath;
    ReadKey key;

    if (graph_->NumReads() == 0) {
        graph_->AddFirstRead(readSequence, &outputPath);
        readPaths_.push_back(outputPath);
        reverseComplemented_.push_back(false);
        key = graph_->NumReads() - 1;
    } else {
        auto c1 = graph_->TryAddRead(readSequence, config, rangeFinder_);
        auto c2 = graph_->TryAddRead(ReverseComplement(readSequence), config, rangeFinder_);

        if (c1->Score() >= c2->Score() && c1->Score() >= minScoreToAdd) {
            graph_->CommitAdd(c1, &outputPath);
            readPaths_.push_back(outputPath);
            reverseComplemented_.push_back(false);
            key = graph_->NumReads() - 1;
        } else if (c2->Score() >= c1->Score() && c2->Score() >= minScoreToAdd) {
            graph_->CommitAdd(c2, &outputPath);
            readPaths_.push_back(outputPath);
            reverseComplemented_.push_back(true);
            key = graph_->NumReads() - 1;
        } else {
            key = -1;
        }

        delete c1;
        delete c2;
    }
    return key;
}

std::shared_ptr<const PoaConsensus> SparsePoa::FindConsensus(
    int minCoverage, std::vector<PoaAlignmentSummary>* summaries) const
{
    AlignConfig config = DefaultPoaConfig(AlignMode::LOCAL);
    std::shared_ptr<const PoaConsensus> pc(graph_->FindConsensus(config, minCoverage));
    std::string css = pc->Sequence;

    if (summaries != NULL) {
        summaries->clear();

        // digest the consensus path consensus into map(vtx, pos)
        // the fold over the readPaths
        std::map<Vertex, size_t> cssPosition;

        int i = 0;
        for (Vertex v : pc->Path) {
            cssPosition[v] = i;
            i++;
        }

        for (size_t readId = 0; readId < graph_->NumReads(); readId++) {
            size_t readS = 0, readE = 0;
            size_t cssS = 0, cssE = 0;
            bool foundStart = false;
            size_t nErr = 0;

            const std::vector<Vertex>& readPath = readPaths_[readId];

            for (size_t readPos = 0; readPos < readPath.size(); readPos++) {
                Vertex v = readPath[readPos];
                if (cssPosition.find(v) != cssPosition.end()) {
                    if (!foundStart) {
                        cssS = cssPosition[v];
                        readS = readPos;
                        foundStart = true;
                    }

                    cssE = cssPosition[v] + 1;
                    readE = readPos + 1;
                } else {
                    nErr += 1;
                }
            }

            Interval readExtent(readS, readE);
            Interval cssExtent(cssS, cssE);

            PoaAlignmentSummary summary;
            summary.ReverseComplementedRead = reverseComplemented_[readId];
            summary.ExtentOnRead = readExtent;
            summary.ExtentOnConsensus = cssExtent;
            summary.AlignmentIdentity = std::max(0.0f, 1.0f - 1.0f * nErr / cssPosition.size());

            (*summaries).push_back(summary);
        }
    }

    return pc;
}

std::string SparsePoa::ToGraphViz(int flags, const PoaConsensus* pc) const
{
    return graph_->ToGraphViz(flags, pc);
}

void SparsePoa::WriteGraphVizFile(const std::string& filename, int flags,
                                  const PoaConsensus* pc) const
{
    graph_->WriteGraphVizFile(filename, flags, pc);
}

void SparsePoa::WriteGraphCsvFile(const std::string& filename) const
{
    graph_->WriteGraphCsvFile(filename);
}

void SparsePoa::PruneGraph(const int minCoverage) { graph_->PruneGraph(minCoverage); }

void SparsePoa::repCheck()
{
    assert(graph_->NumReads() == readPaths_.size());
    assert(graph_->NumReads() == reverseComplemented_.size());
}

}  // namespace Poa
}  // namespace PacBio
