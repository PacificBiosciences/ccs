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

#include <string>
#include <utility>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>

#include <pacbio/consensus/align/AlignConfig.h>
#include <pacbio/consensus/poa/PoaConsensus.h>

namespace PacBio {
namespace Consensus {

AlignConfig DefaultPoaConfig(AlignMode mode)
{
    AlignParams params(3, -5, -4, -4);
    AlignConfig config(params, mode);
    return config;
}

PoaConsensus::PoaConsensus(const std::string& css, const PoaGraph& g,
                           const std::vector<size_t>& cssPath)
    : Sequence(css), Graph(g), Path(cssPath)
{
}

PoaConsensus::PoaConsensus(const std::string& css, const detail::PoaGraphImpl& gi,
                           const std::vector<size_t>& cssPath)
    : Sequence(css), Graph(gi), Path(cssPath)
{
}

PoaConsensus::~PoaConsensus() {}
const PoaConsensus* PoaConsensus::FindConsensus(const std::vector<std::string>& reads)
{
    return FindConsensus(reads, DefaultPoaConfig(AlignMode::GLOBAL), -INT_MAX);
}

const PoaConsensus* PoaConsensus::FindConsensus(const std::vector<std::string>& reads,
                                                const AlignConfig& config, int minCoverage)
{
    PoaGraph pg;
    for (const std::string& read : reads) {
        if (read.length() == 0) {
            throw std::invalid_argument("input sequences must have nonzero length.");
        }
        pg.AddRead(read, config);
    }
    return pg.FindConsensus(config, minCoverage);
}

const PoaConsensus* PoaConsensus::FindConsensus(const std::vector<std::string>& reads,
                                                AlignMode mode, int minCoverage)
{
    return FindConsensus(reads, DefaultPoaConfig(mode), minCoverage);
}

std::string PoaConsensus::ToGraphViz(int flags) const { return Graph.ToGraphViz(flags, this); }
void PoaConsensus::WriteGraphVizFile(std::string filename, int flags) const
{
    Graph.WriteGraphVizFile(filename, flags, this);
}

}  // namespace Consensus
}  // namespace PacBio
