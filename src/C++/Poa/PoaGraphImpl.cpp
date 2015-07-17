#include "PoaGraphImpl.hpp"

#include <ConsensusCore/Align/AlignConfig.hpp>
#include <ConsensusCore/Interval.hpp>
#include <ConsensusCore/Poa/PoaConsensus.hpp>
#include <ConsensusCore/Poa/PoaGraph.hpp>
#include <ConsensusCore/Poa/RangeFinder.hpp>
#include <ConsensusCore/Utils.hpp>

#include <boost/graph/copy.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/graphviz.hpp>

#include <set>
#include <iostream>


namespace boost
{
    using ConsensusCore::detail::VertexInfoMap;
    using ConsensusCore::PoaConsensus;
    using ConsensusCore::PoaGraph;
    using ConsensusCore::detail::VD;
    using boost::format;

    class my_label_writer
    {
    public:
        my_label_writer(VertexInfoMap map, bool color, bool verbose, const PoaConsensus* pc = NULL)
            : map_(map),
              cssVtxs_(),
              color_(color),
              verbose_(verbose)
        {
            if (pc != NULL)
            {
                cssVtxs_.insert(pc->Path.begin(), pc->Path.end());
            }
        }

        template <class descriptor>
        void operator()(std::ostream& out, const descriptor& v) const
        {
            PoaGraph::Vertex vertexId = map_[v].Id;

            std::string nodeColoringAttribute =
                (color_ && isInConsensus(vertexId) ?
                 " style=\"filled\", fillcolor=\"lightblue\" ," : "");

            if (!verbose_)
            {
                out << format("[shape=Mrecord,%s label=\"{ %c | %d }\"]")
                    % nodeColoringAttribute
                    % map_[v].Base
                    % map_[v].Reads;
            }
            else
            {
                out <<  format("[shape=Mrecord,%s label=\"{ "
                               "{ %d | %c } |"
                               "{ %d | %d } |"
                               "{ %0.2f | %0.2f } }\"]")
                    % nodeColoringAttribute
                    % vertexId % map_[v].Base
                    % map_[v].Reads % map_[v].SpanningReads
                    % map_[v].Score % map_[v].ReachingScore;
            }
        }
    private:

        bool isInConsensus(PoaGraph::Vertex v) const
        {
            return cssVtxs_.find(v) != cssVtxs_.end();
        }

        VertexInfoMap map_;
        std::set<PoaGraph::Vertex> cssVtxs_;
        bool color_;
        bool verbose_;
    };
}

namespace ConsensusCore {
namespace detail {

    // ----------------- PoaAlignmentMatrixImpl ---------------------


    PoaAlignmentMatrixImpl::~PoaAlignmentMatrixImpl()
    {
        foreach (AlignmentColumnMap::value_type& kv, columns_)
        {
            delete kv.second;
        }
    }


    float PoaAlignmentMatrixImpl::Score() const
    {
        return score_;
    }


    // ----------------- PoaGraphImpl ---------------------

    PoaGraphImpl::PoaGraphImpl()
        : g_(),
          vertexInfoMap_(get(vertex_info, g_)),
          numReads_(0),
          totalVertices_(0),
          liveVertices_(0)
    {
        enterVertex_ = addVertex('^', 0);
        exitVertex_  = addVertex('$', 0);
    }

    PoaGraphImpl::PoaGraphImpl(const PoaGraphImpl& other)
        : g_(other.g_),
          vertexInfoMap_(get(vertex_info, g_)),
          enterVertex_(other.enterVertex_),
          exitVertex_(other.exitVertex_),
          numReads_(other.numReads_)
    {}

    PoaGraphImpl::~PoaGraphImpl()
    {}

    void PoaGraphImpl::repCheck() const
    {
        // assert the representation invariant for the object
        foreach (VD v, vertices(g_))
        {
            if (v == enterVertex_)
            {
                assert(in_degree(v, g_) == 0);
                assert(out_degree(v, g_) > 0 || NumReads() == 0);
            }
            else if (v == exitVertex_)
            {
                assert(in_degree(v, g_) > 0 || NumReads() == 0);
                assert(out_degree(v, g_) == 0);
            }
            else
            {
                assert(in_degree(v, g_) > 0);
                assert(out_degree(v, g_) > 0);
            }
        }
    }

    static inline vector<const AlignmentColumn*>
    getPredecessorColumns(const BoostGraph& g,
                          VD v,
                          const AlignmentColumnMap& colMap)
    {
        vector<const AlignmentColumn*> predecessorColumns;
        const AlignmentColumn* predCol;
        foreach (ED e, inEdges(v, g))
        {
            VD u = source(e, g);
            predCol = colMap.at(u);
            assert(predCol != NULL);
            predecessorColumns.push_back(predCol);
        }
        return predecessorColumns;
    }

    PoaConsensus*
    PoaGraphImpl::FindConsensus(const AlignConfig& config, int minCoverage)
    {
        std::vector<VD> bestPath = consensusPath(config.Mode, minCoverage);
        std::string consensusSequence = sequenceAlongPath(g_, vertexInfoMap_, bestPath);
        PoaConsensus* pc = new PoaConsensus(consensusSequence, *this, externalizePath(bestPath));
        return pc;
    }

    const AlignmentColumn*
    PoaGraphImpl::makeAlignmentColumnForExit(VD v,
                                             const AlignmentColumnMap& colMap,
                                             const std::string& sequence,
                                             const AlignConfig& config) const
    {
        assert(out_degree(v, g_) == 0);

        // this is kind of unnecessary as we are only actually using one entry in this column
        int I = sequence.length();
        AlignmentColumn* curCol = new AlignmentColumn(v, I + 1);

        float bestScore = -FLT_MAX;
        VD prevVertex = null_vertex;

        // Under local or semiglobal alignment the vertex $ can be
        // "reached" in the dynamic programming from any other vertex
        // in one step via the End move--not just its predecessors in
        // the graph.  In local alignment, it may have been from any
        // row, not necessarily I.
        if (config.Mode == SEMIGLOBAL || config.Mode == LOCAL)
        {
            foreach (VD u, vertices(g_))
            {
                if (u != exitVertex_)
                {
                    const AlignmentColumn* predCol = colMap.at(u);
                    int prevRow = (config.Mode == LOCAL ? ArgMax(predCol->Score) : I);

                    if (predCol->Score[prevRow] > bestScore)
                    {
                        bestScore = predCol->Score[prevRow];
                        prevVertex = predCol->CurrentVertex;
                    }
                }
            }
        }
        else
        {
            // regular predecessors
            vector<const AlignmentColumn*> predecessorColumns  =
                    getPredecessorColumns(g_, v, colMap);
            foreach (const AlignmentColumn * predCol, predecessorColumns)
            {
                if (predCol->Score[I] > bestScore)
                {
                    bestScore = predCol->Score[I];
                    prevVertex = predCol->CurrentVertex;
                }
            }
        }
        assert(prevVertex != null_vertex);
        curCol->Score[I] = bestScore;
        curCol->PreviousVertex[I] = prevVertex;
        curCol->ReachingMove[I] = EndMove;
        return curCol;
    }

    const AlignmentColumn*
    PoaGraphImpl::makeAlignmentColumn(VD v,
                                      const AlignmentColumnMap& colMap,
                                      const std::string& sequence,
                                      const AlignConfig& config,
                                      int beginRow,
                                      int endRow) const
    {
        AlignmentColumn* curCol = new AlignmentColumn(v, sequence.length() + 1);
        const PoaNode& vertexInfo = vertexInfoMap_[v];
        vector<const AlignmentColumn*> predecessorColumns =
                getPredecessorColumns(g_, v, colMap);

        //
        // handle row 0 separately:
        //
        if (predecessorColumns.size() == 0)
        {
            // if this vertex doesn't have any in-edges it is ^; has
            // no reaching move
            assert(v == enterVertex_);
            curCol->Score[0] = 0;
            curCol->ReachingMove[0] = InvalidMove;
            curCol->PreviousVertex[0] = null_vertex;
        }
        else if (config.Mode == SEMIGLOBAL  || config.Mode == LOCAL)
        {
            // under semiglobal or local alignment, we use the Start move
            curCol->Score[0] = 0;
            curCol->ReachingMove[0] = StartMove;
            curCol->PreviousVertex[0] = enterVertex_;
        }
        else
        {
            // otherwise it's a deletion
            float candidateScore;
            float bestScore = -FLT_MAX;
            VD prevVertex = null_vertex;
            MoveType reachingMove = InvalidMove;

            foreach (const AlignmentColumn * prevCol, predecessorColumns)
            {
                candidateScore = prevCol->Score[0] + config.Params.Delete;
                if (candidateScore > bestScore)
                {
                    bestScore = candidateScore;
                    prevVertex = prevCol->CurrentVertex;
                    reachingMove = DeleteMove;
                }
            }
            assert(reachingMove != InvalidMove);
            curCol->Score[0] = bestScore;
            curCol->ReachingMove[0] = reachingMove;
            curCol->PreviousVertex[0] = prevVertex;
        }

        //
        // tackle remainder of read.
        //
        // i represents position in array
        // readPos=i-1 represents position in read
        for (unsigned int i = 1, readPos = 0;  i <= sequence.length(); i++, readPos++)
        {
            float candidateScore, bestScore;
            VD prevVertex;
            MoveType reachingMove;

            if (config.Mode == LOCAL)
            {
                bestScore = 0;
                prevVertex = enterVertex_;
                reachingMove = StartMove;
            }
            else
            {
                bestScore = -FLT_MAX;
                prevVertex = null_vertex;
                reachingMove = InvalidMove;
            }

            foreach (const AlignmentColumn* prevCol, predecessorColumns)
            {
                // Incorporate (Match or Mismatch)
                bool isMatch = sequence[readPos] == vertexInfo.Base;
                candidateScore = prevCol->Score[i - 1] + (isMatch ?
                                                             config.Params.Match :
                                                             config.Params.Mismatch);
                if (candidateScore > bestScore)
                {
                    bestScore = candidateScore;
                    prevVertex = prevCol->CurrentVertex;
                    reachingMove = (isMatch ? MatchMove : MismatchMove);
                }
                // Delete
                candidateScore = prevCol->Score[i] + config.Params.Delete;
                if (candidateScore > bestScore)
                {
                    bestScore = candidateScore;
                    prevVertex = prevCol->CurrentVertex;
                    reachingMove = DeleteMove;
                }
            }
            // Extra
            candidateScore = curCol->Score[i - 1] + config.Params.Insert;
            if (candidateScore > bestScore)
            {
                bestScore = candidateScore;
                prevVertex = v;
                reachingMove = ExtraMove;
            }
            assert(reachingMove != InvalidMove);
            curCol->Score[i] = bestScore;
            curCol->ReachingMove[i] = reachingMove;
            curCol->PreviousVertex[i] = prevVertex;
        }

        return curCol;
    }

    void PoaGraphImpl::AddRead(const std::string& readSeq,
                               const AlignConfig& config,
                               SdpRangeFinder* rangeFinder,
                               std::vector<Vertex>* readPathOutput)
    {
        if (NumReads() == 0)
        {
            AddFirstRead(readSeq, readPathOutput);
        }
        else
        {
            PoaAlignmentMatrixImpl* mat = TryAddRead(readSeq, config, rangeFinder);
            CommitAdd(mat, readPathOutput);
            delete mat;
        }
    }

    void PoaGraphImpl::AddFirstRead(const std::string& readSeq,
                                    std::vector<Vertex>* readPathOutput)
    {
        DEBUG_ONLY(repCheck());
        assert(readSeq.length() > 0);
        assert(numReads_ == 0);

        threadFirstRead(readSeq, readPathOutput);
        numReads_++;

        DEBUG_ONLY(repCheck());
    }

    PoaAlignmentMatrixImpl*
    PoaGraphImpl::TryAddRead(const std::string& readSeq,
                             const AlignConfig& config,
                             SdpRangeFinder* rangeFinder) const
    {
        DEBUG_ONLY(repCheck());
        assert(readSeq.length() > 0);
        assert(numReads_ > 0);

        // Prepare the range finder, if applicable
        if (rangeFinder != NULL)
        {
            // NB: no minCoverage applicable here; this
            // "intermediate" consensus may include extra sequence
            // at either end
            std::vector<VD> cssPath = consensusPath(config.Mode);
            std::string cssSeq = sequenceAlongPath(g_, vertexInfoMap_, cssPath);
            rangeFinder->InitRangeFinder(*this, externalizePath(cssPath), cssSeq, readSeq);
        }

        // Calculate alignment columns of sequence vs. graph, using sparsity if
        // we have a range finder.
        PoaAlignmentMatrixImpl* mat = new PoaAlignmentMatrixImpl();
        mat->readSequence_ = readSeq;
        mat->mode_ = config.Mode;

        vector<VD> sortedVertices(num_vertices(g_));
        topological_sort(g_, sortedVertices.rbegin());
        const AlignmentColumn* curCol;
        foreach (VD v, sortedVertices)
        {
            if (v != exitVertex_)
            {
                Interval rowRange;
                if (rangeFinder) {
                    rowRange = rangeFinder->FindAlignableRange(externalize(v));
                } else {
                    rowRange = Interval(0, readSeq.size());
                }
                curCol = makeAlignmentColumn(v, mat->columns_, readSeq, config, rowRange.Begin, rowRange.End);
            }
            else {
                curCol = makeAlignmentColumnForExit(v, mat->columns_, readSeq, config);
            }
            mat->columns_[v] = curCol;
        }

        mat->score_ = mat->columns_[exitVertex_]->Score[readSeq.size()];
        DEBUG_ONLY(repCheck());

        return mat;
    }

    void
    PoaGraphImpl::CommitAdd(PoaAlignmentMatrix* mat_, std::vector<Vertex>* readPathOutput)
    {
        DEBUG_ONLY(repCheck());

        PoaAlignmentMatrixImpl* mat = static_cast<PoaAlignmentMatrixImpl*>(mat_);
        tracebackAndThread(mat->readSequence_, mat->columns_, mat->mode_, readPathOutput);
        numReads_++;

        DEBUG_ONLY(repCheck());
    }

    size_t PoaGraphImpl::NumReads() const
    {
       return numReads_;
    }

    string PoaGraphImpl::ToGraphViz(int flags, const PoaConsensus* pc) const
    {
       std::stringstream ss;
       write_graphviz(ss, g_, my_label_writer(vertexInfoMap_,
                                              flags & PoaGraph::COLOR_NODES,
                                              flags & PoaGraph::VERBOSE_NODES,
                                              pc));
       return ss.str();
    }

   void
   PoaGraphImpl::WriteGraphVizFile(string filename, int flags, const PoaConsensus* pc) const
   {
       std::ofstream outfile(filename.c_str());
       outfile << ToGraphViz(flags, pc);
       outfile.close();
   }

}} // ConsensusCore::detail
