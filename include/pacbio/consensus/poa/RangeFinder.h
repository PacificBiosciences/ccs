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

#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <pacbio/consensus/poa/PoaGraph.h>

namespace PacBio {
namespace Consensus {
namespace detail {

// an Anchor represents a point (cssPos, readPos)
typedef std::pair<size_t, size_t> SdpAnchor;
typedef std::vector<SdpAnchor> SdpAnchorVector;

class PoaGraphImpl;

// using PacBio::Consensus::PoaGraph;

//
// SdpRangeFinder objects are responsible for identifying the range
// of read positions that we should seek to align to a POA vertex;
// this implementation uses SDP to identify fairly narrow bands,
// enabling sparse memory usage.
//
// This is an abstract class that will be inherited in a client
// library that has access to an SDP method.
//
// RangeFinder state goes away on next call to InitRangeFinder.  We could
// have dealt with this using a factory pattern but bleh.
//
class SdpRangeFinder
{
private:
    std::map<PoaGraph::Vertex, std::pair<int, int>> alignableReadIntervalByVertex_;

public:
    virtual ~SdpRangeFinder();

    void InitRangeFinder(const PoaGraphImpl& poaGraph,
                         const std::vector<PoaGraph::Vertex>& consensusPath,
                         const std::string& consensusSequence, const std::string& readSequence);

    // TODO: write contract
    std::pair<int, int> FindAlignableRange(PoaGraph::Vertex v);

protected:
    // TODO: write contract
    virtual SdpAnchorVector FindAnchors(const std::string& consensusSequence,
                                        const std::string& readSequence) const = 0;
};

}  // namespace detail
}  // namespace Consensus
}  // namespace PacBio
