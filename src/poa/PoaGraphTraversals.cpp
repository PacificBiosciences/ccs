// Copyright (c) 2011-2015, Pacific Biosciences of California, Inc.
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

#include <boost/foreach.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/unordered_set.hpp>

#include <pacbio/consensus/poa/PoaGraph.h>

#include "PoaGraphImpl.h"
#include "VectorL.h"

namespace PacBio {
namespace Consensus {
namespace detail {

std::string sequenceAlongPath(const BoostGraph& g, const VertexInfoMap& vertexInfoMap,
                              const std::vector<VD>& path)
{
    std::stringstream ss;
    for (const VD v : path) {
        ss << vertexInfoMap[v].Base;
    }
    return ss.str();
}

boost::unordered_set<VD> SpanningDFS(const VD start, const VD end, const BoostGraph& g)
{
    std::vector<VD> stack;
    boost::unordered_set<VD> fwd;
    boost::unordered_set<VD> rev;
    // find all vertices reachable from start
    stack.push_back(start);
    do {
        VD v = stack.back();
        stack.pop_back();
        // mark those we've already visited,
        //   if so, skip
        if (fwd.find(v) != fwd.end()) continue;
        fwd.insert(v);
        BOOST_FOREACH (const ED& e, out_edges(v, g)) {
            stack.push_back(target(e, g));
        }
    } while (!stack.empty());
    // find all vertices that can reach end
    //   FROM start
    stack.push_back(end);
    do {
        VD v = stack.back();
        stack.pop_back();
        // if it's not been visited in the forward pass,
        //   or we've already visited it here, skip
        if (fwd.find(v) == fwd.end() || rev.find(v) != rev.end()) continue;
        rev.insert(v);
        BOOST_FOREACH (const ED& e, in_edges(v, g)) {
            stack.push_back(source(e, g));
        }
    } while (!stack.empty());
    return rev;
}

std::vector<VD> PoaGraphImpl::sortedVertices() const
{
    vector<VD> sv(num_vertices(g_));
    topological_sort(g_, sv.rbegin());
    return sv;
}

void PoaGraphImpl::tagSpan(VD start, VD end)
{
    boost::unordered_set<VD> vertices = SpanningDFS(start, end, g_);
    for (const VD v : vertices) {
        vertexInfoMap_[v].SpanningReads++;
    }
}

std::vector<VD> PoaGraphImpl::consensusPath(AlignMode mode, int minCoverage) const
{
    // Pat's note on the approach here:
    //
    // "A node gets a score of NumReads if all reads go through
    //  it, and a score of -NumReads if no reads go through it The
    //  shift of -0.0001 breaks ties in favor of skipping
    //  half-full nodes.  In the 2 reads case this will get rid of
    //  insertions which are the more common error."
    //
    // The interpretation of minCoverage (which is applicable only
    // for LOCAL, SEMIGLOBAL modes) is that it represents
    // application-specific knowledge of the basal coverage level
    // of reads in the template, such that if a node is contained
    // in fewer than minCoverage reads, it will be penalized
    // against inclusion in the consensus.
    int totalReads = NumReads();

    std::list<VD> path;
    std::list<VD> sortedVertices(num_vertices(g_));
    topological_sort(g_, sortedVertices.rbegin());
    unordered_map<VD, VD> bestPrevVertex;

    // ignore ^ and $
    // TODO(dalexander): find a cleaner way to do this
    vertexInfoMap_[sortedVertices.front()].ReachingScore = 0;
    sortedVertices.pop_back();
    sortedVertices.pop_front();

    VD bestVertex = null_vertex;
    float bestReachingScore = -FLT_MAX;
    for (const VD v : sortedVertices) {
        PoaNode& vInfo = vertexInfoMap_[v];
        int containingReads = vInfo.Reads;
        int spanningReads = vInfo.SpanningReads;
        float score =
            (mode != AlignMode::GLOBAL)
                ? (2 * containingReads - 1 * std::max(spanningReads, minCoverage) - 0.0001f)
                : (2 * containingReads - 1 * totalReads - 0.0001f);
        vInfo.Score = score;
        vInfo.ReachingScore = score;
        bestPrevVertex[v] = null_vertex;
        for (const ED& e : inEdges(v, g_)) {
            VD sourceVertex = source(e, g_);
            float rsc = score + vertexInfoMap_[sourceVertex].ReachingScore;
            if (rsc > vInfo.ReachingScore) {
                vInfo.ReachingScore = rsc;
                bestPrevVertex[v] = sourceVertex;
            }
            if (rsc > bestReachingScore) {
                bestVertex = v;
                bestReachingScore = rsc;
            }
            // if the score is the same, the order we've encountered vertices
            //   might not be deterministic. Fix this by comparing on
            //   vertex_index
            else if (rsc == bestReachingScore) {
                if (get(vertex_index, g_, v) < get(vertex_index, g_, bestVertex)) bestVertex = v;
            }
        }
    }
    assert(bestVertex != null_vertex);

    // trace back from best-scoring vertex
    VD v = bestVertex;
    while (v != null_vertex) {
        path.push_front(v);
        v = bestPrevVertex[v];
    }
    return std::vector<VD>(path.begin(), path.end());
}

void PoaGraphImpl::threadFirstRead(std::string sequence, std::vector<Vertex>* outputPath)
{
    // first sequence in the alignment
    VD u = null_vertex, v;
    VD startSpanVertex = null_vertex, endSpanVertex;
    int readPos = 0;

    if (outputPath) {
        outputPath->clear();
    }

    for (const char base : sequence) {
        v = addVertex(base);
        if (outputPath) {
            outputPath->push_back(externalize(v));
        }
        if (readPos == 0) {
            add_edge(enterVertex_, v, g_);
            startSpanVertex = v;
        } else {
            add_edge(u, v, g_);
        }
        u = v;
        readPos++;
    }
    assert(startSpanVertex != null_vertex);
    assert(u != null_vertex);
    endSpanVertex = u;
    add_edge(u, exitVertex_, g_);  // terminus -> $
    tagSpan(startSpanVertex, endSpanVertex);
}

void PoaGraphImpl::tracebackAndThread(std::string sequence,
                                      const AlignmentColumnMap& alignmentColumnForVertex,
                                      AlignMode alignMode, std::vector<Vertex>* outputPath)
{
    const int I = sequence.length();

    // perform traceback from (I,$), threading the new sequence into
    // the graph as we go.
    int i = I;
    const AlignmentColumn* curCol;
    VD v = null_vertex, forkVertex = null_vertex;
    VD u = exitVertex_;
    VD startSpanVertex;
    VD endSpanVertex = alignmentColumnForVertex.at(exitVertex_)->PreviousVertex[I];

    if (outputPath) {
        outputPath->resize(I);
        std::fill(outputPath->begin(), outputPath->end(), (size_t)-1);
    }

#define READPOS (i - 1)
#define VERTEX_ON_PATH(readPos, v)                 \
    if (outputPath) {                              \
        (*outputPath)[(readPos)] = externalize(v); \
    }

    while (!(u == enterVertex_ && i == 0)) {
        // u -> v
        // u: current vertex
        // v: vertex last visited in traceback (could be == u)
        // forkVertex: the vertex that will be the target of a new edge

        Vertex uExt = this->externalize(u);  // DEBUGGING
        Vertex vExt = this->externalize(v);  // DEBUGGING

        curCol = alignmentColumnForVertex.at(u);
        assert(curCol != NULL);
        PoaNode& curNodeInfo = vertexInfoMap_[u];
        VD prevVertex = curCol->PreviousVertex[i];
        MoveType reachingMove = curCol->ReachingMove[i];

        if (reachingMove == StartMove) {
            assert(v != null_vertex);

            if (forkVertex == null_vertex) {
                forkVertex = v;
            }
            // In local model thread read bases, adjusting i (should stop at 0)
            while (i > 0) {
                assert(alignMode == AlignMode::LOCAL);
                VD newForkVertex = addVertex(sequence[READPOS]);
                add_edge(newForkVertex, forkVertex, g_);
                VERTEX_ON_PATH(READPOS, newForkVertex);
                forkVertex = newForkVertex;
                i--;
            }
        } else if (reachingMove == EndMove) {
            assert(forkVertex == null_vertex && u == exitVertex_ && v == null_vertex);

            forkVertex = exitVertex_;

            if (alignMode == AlignMode::LOCAL) {
                // Find the row # we are coming from, walk
                // back to there, threading read bases onto
                // graph via forkVertex, adjusting i.
                const AlignmentColumn* prevCol = alignmentColumnForVertex.at(prevVertex);
                int prevRow = ArgMax(prevCol->Score);

                while (i > static_cast<int>(prevRow)) {
                    VD newForkVertex = addVertex(sequence[READPOS]);
                    add_edge(newForkVertex, forkVertex, g_);
                    VERTEX_ON_PATH(READPOS, newForkVertex);
                    forkVertex = newForkVertex;
                    i--;
                }
            }
        } else if (reachingMove == MatchMove) {
            VERTEX_ON_PATH(READPOS, u);
            // if there is an extant forkVertex, join it
            if (forkVertex != null_vertex) {
                add_edge(u, forkVertex, g_);
                forkVertex = null_vertex;
            }
            // add to existing node
            curNodeInfo.Reads++;
            i--;
        } else if (reachingMove == DeleteMove) {
            if (forkVertex == null_vertex) {
                forkVertex = v;
            }
        } else if (reachingMove == ExtraMove || reachingMove == MismatchMove) {
            // begin a new arc with this read base
            VD newForkVertex = addVertex(sequence[READPOS]);
            if (forkVertex == null_vertex) {
                forkVertex = v;
            }
            add_edge(newForkVertex, forkVertex, g_);
            VERTEX_ON_PATH(READPOS, newForkVertex);
            forkVertex = newForkVertex;
            i--;
        } else {
            throw std::runtime_error("unreachable");
        }

        v = u;
        u = prevVertex;
    }
    startSpanVertex = v;

    // if there is an extant forkVertex, join it to enterVertex
    if (forkVertex != null_vertex) {
        add_edge(enterVertex_, forkVertex, g_);
        startSpanVertex = forkVertex;
        forkVertex = null_vertex;
    }

    if (startSpanVertex != exitVertex_) {
        tagSpan(startSpanVertex, endSpanVertex);
    }

    // all filled in?
    assert(outputPath == NULL ||
           std::find(outputPath->begin(), outputPath->end(), ((size_t)-1)) == outputPath->end());

#undef READPOS
#undef VERTEX_ON_PATH
}

static boost::unordered_set<VD> childVertices(VD v, const BoostGraph& g)
{
    boost::unordered_set<VD> result;
    BOOST_FOREACH (const ED& e, out_edges(v, g)) {
        result.insert(target(e, g));
    }
    return result;
}

static boost::unordered_set<VD> parentVertices(VD v, const BoostGraph& g)
{
    boost::unordered_set<VD> result;
    BOOST_FOREACH (const ED& e, in_edges(v, g)) {
        result.insert(source(e, g));
    }
    return result;
}

vector<ScoredMutation>* PoaGraphImpl::findPossibleVariants(
    const std::vector<Vertex>& bestPath) const
{
    std::vector<VD> bestPath_ = internalizePath(bestPath);

    // Return value will be deallocated by PoaConsensus destructor.
    vector<ScoredMutation>* variants = new vector<ScoredMutation>();

    for (int i = 2; i < (int)bestPath_.size() - 2; i++)  // NOLINT
    {
        VD v = bestPath_[i];
        boost::unordered_set<VD> children = childVertices(v, g_);

        // Look for a direct edge from the current node to the node
        // two spaces down---suggesting a deletion with respect to
        // the consensus sequence.
        if (children.find(bestPath_[i + 2]) != children.end()) {
            float score = -vertexInfoMap_[bestPath_[i + 1]].Score;
            variants->push_back(Mutation(MutationType::DELETION, i + 1, '-').WithScore(score));
        }

        // Look for a child node that connects immediately back to i + 1.
        // This indicates we should try inserting the base at i + 1.

        // Parents of (i + 1)
        boost::unordered_set<VD> lookBack = parentVertices(bestPath_[i + 1], g_);

        // (We could do this in STL using std::set sorted on score, which would
        // then
        // provide an intersection mechanism (in <algorithm>) but that actually
        // ends
        // up being more code.  Sad.)
        float bestInsertScore = -FLT_MAX;
        VD bestInsertVertex = null_vertex;

        for (const VD v : children) {
            boost::unordered_set<VD>::iterator found = lookBack.find(v);
            if (found != lookBack.end()) {
                float score = vertexInfoMap_[*found].Score;
                if (score > bestInsertScore) {
                    bestInsertScore = score;
                    bestInsertVertex = *found;
                } else if (score == bestInsertScore) {
                    if (get(vertex_index, g_, *found) < get(vertex_index, g_, bestInsertVertex))
                        bestInsertVertex = *found;
                }
            }
        }

        if (bestInsertVertex != null_vertex) {
            char base = vertexInfoMap_[bestInsertVertex].Base;
            variants->push_back(
                Mutation(MutationType::INSERTION, i + 1, base).WithScore(bestInsertScore));
        }

        // Look for a child node not in the consensus that connects immediately
        // to i + 2.  This indicates we should try mismatching the base i + 1.

        // Parents of (i + 2)
        lookBack = parentVertices(bestPath_[i + 2], g_);

        float bestMismatchScore = -FLT_MAX;
        VD bestMismatchVertex = null_vertex;

        for (const VD v : children) {
            if (v == bestPath_[i + 1]) continue;

            boost::unordered_set<VD>::iterator found = lookBack.find(v);
            if (found != lookBack.end()) {
                float score = vertexInfoMap_[*found].Score;
                if (score > bestMismatchScore) {
                    bestMismatchScore = score;
                    bestMismatchVertex = *found;
                } else if (score == bestMismatchScore) {
                    if (get(vertex_index, g_, *found) < get(vertex_index, g_, bestMismatchVertex))
                        bestMismatchVertex = *found;
                }
            }
        }

        if (bestMismatchVertex != null_vertex) {
            // TODO(dalexander): As implemented (compatibility), this returns
            // the score of the mismatch node. I think it should return the
            // score
            // difference, no?
            char base = vertexInfoMap_[bestMismatchVertex].Base;
            variants->push_back(
                Mutation(MutationType::SUBSTITUTION, i + 1, base).WithScore(bestMismatchScore));
        }
    }
    return variants;
}

}  // namespace detail
}  // namespace Consensus
}  // namespace PacBio
