// Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.
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

#include <climits>
#include <string>
#include <utility>
#include <vector>

#include <pacbio/consensus/align/AlignConfig.h>

namespace PacBio {
namespace Consensus {

// fwd decls
namespace detail {

class PoaGraphImpl;
class SdpRangeFinder;

}  // namespace detail

struct PoaConsensus;

class PoaAlignmentMatrix
{
public:
    virtual ~PoaAlignmentMatrix(){};
    virtual float Score() const = 0;
    virtual size_t NumRows() const = 0;
    virtual size_t NumCols() const = 0;
    virtual void Print() const = 0;
};

/// \brief An object representing a Poa (partial-order alignment) graph
class PoaGraph
{
public:
    typedef size_t Vertex;
    typedef size_t ReadId;

    static const Vertex NullVertex = (Vertex)-1;

public:  // Flags enums for specifying GraphViz output features
    enum
    {
        COLOR_NODES = 0x1,
        VERBOSE_NODES = 0x2
    };

public:
    PoaGraph();
    PoaGraph(const PoaGraph& other);
    PoaGraph(const detail::PoaGraphImpl& o);  // NB: this performs a copy
    ~PoaGraph();

    //
    // Easy API
    //
    void AddRead(const std::string& sequence, const AlignConfig& config,
                 detail::SdpRangeFinder* rangeFinder = NULL,
                 std::vector<Vertex>* readPathOutput = NULL);

    //
    // API for more control
    //
    void AddFirstRead(const std::string& sequence, std::vector<Vertex>* readPathOutput = NULL);

    PoaAlignmentMatrix* TryAddRead(const std::string& sequence, const AlignConfig& config,
                                   detail::SdpRangeFinder* rangeFinder = NULL) const;

    void CommitAdd(PoaAlignmentMatrix* mat, std::vector<Vertex>* readPathOutput = NULL);

    // ----------

    size_t NumReads() const;

    std::string ToGraphViz(int flags = 0, const PoaConsensus* pc = NULL) const;

    void WriteGraphVizFile(const std::string& filename, int flags = 0,
                           const PoaConsensus* pc = NULL) const;

    const PoaConsensus* FindConsensus(const AlignConfig& config, int minCoverage = -INT_MAX) const;

private:
    detail::PoaGraphImpl* impl;
};

}  // namespace Consensus
}  // namespace PacBio
