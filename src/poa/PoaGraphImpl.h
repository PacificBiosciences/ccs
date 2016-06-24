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

#include <cfloat>
#include <climits>
#include <vector>

#include <boost/config.hpp>
#include <boost/format.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>

#include <pacbio/consensus/Mutation.h>
#include <pacbio/consensus/align/AlignConfig.h>
#include <pacbio/consensus/poa/PoaGraph.h>

#include "BoostGraph.h"
#include "PoaAlignmentMatrix.h"

using std::string;
using std::vector;

///
///  Boost graph library typedefs, properties, and graphviz output.
///
using namespace boost;  // NOLINT

namespace PacBio {
namespace Consensus {
namespace detail {

// FWD
class SdpRangeFinder;

struct PoaNode
{
    PoaGraph::Vertex Id;
    char Base;
    int Reads;
    // move the below out of here?
    int SpanningReads;
    float Score;
    float ReachingScore;

    void Init(size_t id, char base, int reads)
    {
        this->Id = id;
        this->Base = base;
        this->Reads = reads;
        this->SpanningReads = 0;
        this->Score = 0;
        this->ReachingScore = 0;
    }

    PoaNode() { Init(0, 'N', 0); }
    PoaNode(size_t id, char base) { Init(id, base, 1); }
    PoaNode(size_t id, char base, int reads) { Init(id, base, reads); }
};

// External-facing vertex id type
typedef size_t Vertex;

struct EdgeComparator
{
    EdgeComparator(const BoostGraph& g) : g_(g) {}
    // want lex. comparison... just using pair to get it..
    bool operator()(ED e1, ED e2) const
    {
        std::pair<int, int> vt1, vt2;
        vt1 = std::make_pair(get(vertex_index, g_, source(e1, g_)),
                             get(vertex_index, g_, target(e1, g_)));
        vt2 = std::make_pair(get(vertex_index, g_, source(e2, g_)),
                             get(vertex_index, g_, target(e2, g_)));
        return (vt1 < vt2);
    }

private:
    const BoostGraph& g_;
};

inline std::vector<ED> inEdges(VD v, const BoostGraph& g)
{
    // This is a sad workaround the nondeterministic order of iteration
    // from BGL's in_edges. (see: http://stackoverflow.com/questions/30968690/)

    // Unfortunately, we can't just use the boost::sort range adapter
    // because it requires an underlying random access iterator, which
    // we can't get from the std::set container.
    graph_traits<BoostGraph>::in_edge_iterator begin, end;
    tie(begin, end) = in_edges(v, g);
    std::vector<ED> tmp(begin, end);
    std::sort(tmp.begin(), tmp.end(), EdgeComparator(g));
    return tmp;
}

class PoaGraphImpl
{
    friend class SdpRangeFinder;

    BoostGraph g_;
    VertexInfoMap vertexInfoMap_;  // NB: this is a reference type and refers to
                                   // an "internal" property map
    index_map_t indexMap_;
    VD enterVertex_;
    VD exitVertex_;
    size_t numReads_;
    size_t totalVertices_;               // includes "ex"-vertices which have since been removed
    size_t liveVertices_;                // vertices that are in the graph.  this is needed
                                         // for algorithms.
    std::map<Vertex, VD> vertexLookup_;  // external ID -> internal ID

    void repCheck() const;

    VD addVertex(char base, int nReads = 1)
    {
        VD vd = add_vertex(g_);
        Vertex vExt = totalVertices_++;
        vertexInfoMap_[vd] = PoaNode(vExt, base, nReads);
        vertexLookup_[vExt] = vd;
        indexMap_[vd] = liveVertices_++;
        return vd;
    }

    //
    // utility routines
    //
    const AlignmentColumn* makeAlignmentColumn(VD v,
                                               const AlignmentColumnMap& alignmentColumnForVertex,
                                               const std::string& sequence,
                                               const AlignConfig& config, int beginRow,
                                               int endRow) const;

    const AlignmentColumn* makeAlignmentColumnForExit(
        VD v, const AlignmentColumnMap& alignmentColumnForVertex, const std::string& sequence,
        const AlignConfig& config) const;

public:
    //
    // Vertex id translation
    //

    Vertex externalize(VD vd) const
    {
        if (vd == null_vertex) {
            return PoaGraph::NullVertex;
        } else {
            return vertexInfoMap_[vd].Id;
        }
    }

    VD internalize(Vertex vertex) const
    {
        if (vertex == PoaGraph::NullVertex) {
            return null_vertex;
        }
        return vertexLookup_.at(vertex);
    }

    std::vector<Vertex> externalizePath(const std::vector<VD>& vds) const
    {
        std::vector<Vertex> out(vds.size(), 0);
        for (size_t i = 0; i < vds.size(); i++) {
            out[i] = externalize(vds[i]);
        }
        return out;
    }

    std::vector<VD> internalizePath(const std::vector<Vertex>& vertices) const
    {
        std::vector<VD> out(vertices.size(), null_vertex);
        for (size_t i = 0; i < vertices.size(); i++) {
            out[i] = internalize(vertices[i]);
        }
        return out;
    }

public:
    //
    // POA node lookup
    //
    const PoaNode& getPoaNode(VD v) const { return vertexInfoMap_[v]; }
public:
    //
    // Graph traversal functions, defined in PoaGraphTraversals
    //

    std::vector<VD> sortedVertices() const;

    void tagSpan(VD start, VD end);

    std::vector<VD> consensusPath(AlignMode mode, int minCoverage = -INT_MAX) const;

    void threadFirstRead(std::string sequence, std::vector<Vertex>* readPathOutput = NULL);

    void tracebackAndThread(std::string sequence,
                            const AlignmentColumnMap& alignmentColumnForVertex, AlignMode mode,
                            std::vector<Vertex>* readPathOutput = NULL);

    vector<ScoredMutation>* findPossibleVariants(const std::vector<Vertex>& bestPath) const;

public:
    PoaGraphImpl();
    PoaGraphImpl(const PoaGraphImpl& other);
    ~PoaGraphImpl();

    void AddRead(const std::string& sequence, const AlignConfig& config,
                 SdpRangeFinder* rangeFinder = NULL, std::vector<Vertex>* readPathOutput = NULL);

    void AddFirstRead(const std::string& sequence, std::vector<Vertex>* readPathOutput = NULL);

    PoaAlignmentMatrix* TryAddRead(const std::string& sequence, const AlignConfig& config,
                                   SdpRangeFinder* rangeFinder = NULL) const;

    void CommitAdd(PoaAlignmentMatrix* mat, std::vector<Vertex>* readPathOutput = NULL);

    PoaConsensus* FindConsensus(const AlignConfig& config, int minCoverage = -INT_MAX);

    size_t NumReads() const;
    string ToGraphViz(int flags, const PoaConsensus* pc) const;
    void WriteGraphVizFile(const string& filename, int flags, const PoaConsensus* pc) const;
};

// free functions, we should put these all in traversals
std::string sequenceAlongPath(const BoostGraph& g, const VertexInfoMap& vertexInfoMap,
                              const std::vector<VD>& path);

}  // namespace detail
}  // namespace Consensus
}  // namespace PacBio
