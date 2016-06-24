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
// (Based on the original "Partial Order Aligner" by Lee, Grasso, and Sharlow,
//  and an implementation in C# by Patrick Marks)

#include <pacbio/consensus/poa/PoaGraph.h>

#include "PoaGraphImpl.h"

namespace PacBio {
namespace Consensus {

// forward declaration
struct PoaConsensus;

//
// PIMPL idiom delegation
//
void PoaGraph::AddRead(const std::string& sequence, const AlignConfig& config,
                       detail::SdpRangeFinder* rangeFinder, std::vector<Vertex>* readPathOutput)
{
    impl->AddRead(sequence, config, rangeFinder, readPathOutput);
}

void PoaGraph::AddFirstRead(const std::string& sequence, std::vector<Vertex>* readPathOutput)
{
    impl->AddFirstRead(sequence, readPathOutput);
}

PoaAlignmentMatrix* PoaGraph::TryAddRead(const std::string& sequence, const AlignConfig& config,
                                         detail::SdpRangeFinder* rangeFinder) const
{
    return impl->TryAddRead(sequence, config, rangeFinder);
}

void PoaGraph::CommitAdd(PoaAlignmentMatrix* mat, std::vector<Vertex>* readPathOutput)
{
    impl->CommitAdd(mat, readPathOutput);
}

size_t PoaGraph::NumReads() const { return impl->NumReads(); }
const PoaConsensus* PoaGraph::FindConsensus(const AlignConfig& config, int minCoverage) const
{
    return impl->FindConsensus(config, minCoverage);
}

string PoaGraph::ToGraphViz(int flags, const PoaConsensus* pc) const
{
    return impl->ToGraphViz(flags, pc);
}

void PoaGraph::WriteGraphVizFile(const string& filename, int flags, const PoaConsensus* pc) const
{
    impl->WriteGraphVizFile(filename, flags, pc);
}

PoaGraph::PoaGraph() { impl = new detail::PoaGraphImpl(); }
PoaGraph::PoaGraph(const PoaGraph& other) { impl = new detail::PoaGraphImpl(*other.impl); }
PoaGraph::PoaGraph(const detail::PoaGraphImpl& o) { impl = new detail::PoaGraphImpl(o); }
PoaGraph::~PoaGraph() { delete impl; }
}  // namespace Consensus
}  // namespace PacBio
