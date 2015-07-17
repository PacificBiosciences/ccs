
#pragma once
// This header is internal, not part of the API!  ConsensusCore doesn't have a
// good directory-level separation of internal-vs-API, but we can at least prevent
// SWIG export
#ifdef SWIG
#error "PoaGraphImpl.hpp is not an API-facing header!"
#endif  // SWIG

#include <ConsensusCore/Align/AlignConfig.hpp>
#include <ConsensusCore/Poa/PoaGraph.hpp>
#include <ConsensusCore/Matrix/VectorL.hpp>

#include <boost/config.hpp>
#include <boost/format.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/static_assert.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility.hpp>
#include <cfloat>
#include <climits>
#include <vector>

using std::string;
using std::vector;

///
///  Boost graph library typedefs, properties, and graphviz output.
///
using namespace boost; // NOLINT

namespace boost
{
    enum vertex_info_t { vertex_info = 424 };  // a unique #
    BOOST_INSTALL_PROPERTY(vertex, info);
}


namespace ConsensusCore {
namespace detail {

    class SdpRangeFinder;

    enum MoveType
    {
        InvalidMove,  // Invalid move reaching ^ (start)
        StartMove,    // Start move: ^ -> vertex in row 0 of local alignment
        EndMove,      // End move: vertex -> $ in row 0 of local alignment, or
                      //  in global alignment, terminal vertex -> $
        MatchMove,
        MismatchMove,
        DeleteMove,
        ExtraMove
    };

    struct PoaNode
    {
        size_t Id;  // This is the external-facing identifier we use to represent the vertex
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

        PoaNode()
        {
            Init(0, 'N', 0);
        }

        PoaNode(size_t id, char base)
        {
            Init(id, base, 1);
        }

        PoaNode(size_t id, char base, int reads)
        {
            Init(id, base, reads);
        }
    };

    // BGL is intimidating, and it *deserves* your hatred.  But it's
    // the only game in town!
    typedef property<vertex_info_t, PoaNode, property<vertex_index_t, size_t> > vertex_property_t;
    typedef adjacency_list<setS, listS, bidirectionalS, vertex_property_t> BoostGraph;


    // Descriptor types used internally
    typedef graph_traits<BoostGraph>::edge_descriptor   ED;
    typedef graph_traits<BoostGraph>::vertex_descriptor VD;

    // External-facing vertex id type
    typedef size_t Vertex;

    typedef property_map<BoostGraph, vertex_info_t>::type VertexInfoMap;
    typedef property_map<BoostGraph, vertex_index_t>::type index_map_t;
    static const VD null_vertex = graph_traits<BoostGraph>::null_vertex();


    struct EdgeComparator
    {
        EdgeComparator(const BoostGraph& g) : g_(g) {}

        // want lex. comparison... just using pair to get it..
        bool operator() (ED e1, ED e2)
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

    struct AlignmentColumn : noncopyable
    {
        VD CurrentVertex;
        VectorL<float> Score;
        VectorL<MoveType> ReachingMove;
        VectorL<VD> PreviousVertex;

        AlignmentColumn(VD vertex, int len)
            : CurrentVertex(vertex),
              Score(0, len, -FLT_MAX),
              ReachingMove(0, len, InvalidMove),
              PreviousVertex(0, len, null_vertex)
        {}

        AlignmentColumn(VD vertex, int beginRow, int endRow)
            : CurrentVertex(vertex),
              Score(beginRow, endRow, -FLT_MAX),
              ReachingMove(beginRow, endRow, InvalidMove),
              PreviousVertex(beginRow, endRow, null_vertex)
        {}

        ~AlignmentColumn()
        {}

        int BeginRow() const { return Score.BeginRow(); }
        int EndRow()   const { return Score.EndRow();   }
    };


    typedef unordered_map<VD, const AlignmentColumn*> AlignmentColumnMap;

    class PoaAlignmentMatrixImpl : public PoaAlignmentMatrix
    {
    public:
        virtual ~PoaAlignmentMatrixImpl();
        virtual float Score() const;

    public:
        AlignmentColumnMap columns_;
        std::string readSequence_;
        AlignMode mode_;
        float score_;
    };


    class PoaGraphImpl
    {
        friend class SdpRangeFinder;

        BoostGraph g_;
        VertexInfoMap vertexInfoMap_;        // NB: this is a reference type and refers to an "internal" property map
        index_map_t indexMap_;
        VD enterVertex_;
        VD exitVertex_;
        size_t numReads_;
        size_t totalVertices_;               // includes "ex"-vertices which have since been removed
        size_t liveVertices_;                // vertices that are in the graph.  this is needed for algorithms.
        std::map<Vertex, VD> vertexLookup_;  // external ID -> internal ID

        void repCheck() const;

        Vertex externalize(VD vd) const         { return vertexInfoMap_[vd].Id; }
        VD     internalize(Vertex vertex) const { return vertexLookup_.at(vertex); }

        std::vector<Vertex> externalizePath(const std::vector<VD>& vds) const
        {
            std::vector<Vertex> out(vds.size(), 0);
            for (size_t i = 0; i < vds.size(); i++)
            {
                out[i] = externalize(vds[i]);
            }
            return out;
        }

        std::vector<VD> internalizePath(const std::vector<Vertex>& vertices) const
        {
            std::vector<VD> out(vertices.size(), null_vertex);
            for (size_t i = 0; i < vertices.size(); i++)
            {
                out[i] = internalize(vertices[i]);
            }
            return out;
        }

        VD addVertex(char base, int nReads=1)
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
        const AlignmentColumn*
        makeAlignmentColumn(VD v,
                            const AlignmentColumnMap& alignmentColumnForVertex,
                            const std::string& sequence,
                            const AlignConfig& config,
                            int beginRow, int endRow) const;

        const AlignmentColumn*
        makeAlignmentColumnForExit(VD v,
                                   const AlignmentColumnMap& alignmentColumnForVertex,
                                   const std::string& sequence,
                                   const AlignConfig& config) const;

    public:
        //
        // Graph traversal functions, defined in PoaGraphTraversals
        //
        void tagSpan(VD start, VD end);

        std::vector<VD> consensusPath(AlignMode mode, int minCoverage=-INT_MAX) const;

        void threadFirstRead(std::string sequence, std::vector<Vertex>* readPathOutput=NULL);

        void tracebackAndThread
          (std::string sequence,
           const AlignmentColumnMap& alignmentColumnForVertex,
           AlignMode mode,
           std::vector<Vertex>* readPathOutput=NULL);

        vector<ScoredMutation>*
        findPossibleVariants(const std::vector<Vertex>& bestPath) const;

    public:
        PoaGraphImpl();
        PoaGraphImpl(const PoaGraphImpl& other);
        ~PoaGraphImpl();

        void AddRead(const std::string& sequence,
                         const AlignConfig& config,
                         SdpRangeFinder* rangeFinder=NULL,
                         std::vector<Vertex>* readPathOutput=NULL);

        void AddFirstRead(const std::string& sequence,
                              std::vector<Vertex>* readPathOutput=NULL);

        PoaAlignmentMatrixImpl* TryAddRead(const std::string& sequence,
                                           const AlignConfig& config,
                                           SdpRangeFinder* rangeFinder=NULL) const;

        void CommitAdd(PoaAlignmentMatrix* mat, std::vector<Vertex>* readPathOutput=NULL);

        PoaConsensus* FindConsensus(const AlignConfig& config, int minCoverage=-INT_MAX);

        size_t NumReads() const;
        string ToGraphViz(int flags, const PoaConsensus* pc) const;
        void WriteGraphVizFile(string filename, int flags, const PoaConsensus* pc) const;
    };

    // free functions, we should put these all in traversals
    std::string sequenceAlongPath(const BoostGraph& g,
                                  const VertexInfoMap& vertexInfoMap,
                                  const std::vector<VD>& path);

}} // ConsensusCore::detail
