// Author: David Alexander

// The POA "alignment matrix" is a set of alignment columns
// corresponding to each vertex in the graph that the read was aligned
// against.

#pragma once

#include <pacbio/denovo/PoaGraph.h>

#include <boost/unordered_map.hpp>
#include <boost/utility.hpp>
#include <cfloat>

#include "BoostGraph.h"
#include "VectorL.h"

using boost::noncopyable;
using boost::unordered_map;

namespace PacBio {
namespace Poa {
namespace detail {

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

struct AlignmentColumn : noncopyable
{
    VD CurrentVertex;
    VectorL<float> Score;
    VectorL<MoveType> ReachingMove;
    VectorL<VD> PreviousVertex;

    AlignmentColumn(VD vertex, int beginRow, int endRow)
        : CurrentVertex(vertex)
        , Score(beginRow, endRow, -FLT_MAX)
        , ReachingMove(beginRow, endRow, InvalidMove)
        , PreviousVertex(beginRow, endRow, null_vertex)
    {
    }

    ~AlignmentColumn() {}
    size_t BeginRow() const { return Score.BeginRow(); }
    size_t EndRow() const { return Score.EndRow(); }
    bool HasRow(size_t i) const { return (BeginRow() <= i) && (i < EndRow()); }
};

using AlignmentColumnMap = unordered_map<VD, const AlignmentColumn*>;

class PoaAlignmentMatrixImpl : public PoaAlignmentMatrix
{
public:
    virtual ~PoaAlignmentMatrixImpl()
    {
        for (auto& kv : columns_) {
            delete kv.second;
        }
    }

    virtual float Score() const { return score_; }
    size_t NumRows() const { return readSequence_.length() + 1; }
    size_t NumCols() const { return columns_.size(); }
    void Print() const;

public:
    // TODO: why did I leave these public?  why is there no
    // constructor?
    const PoaGraphImpl* graph_;
    AlignmentColumnMap columns_;
    std::string readSequence_;
    PacBio::Align::AlignMode mode_;
    float score_;
};
}
}
}  // PacBio::Poa::detail
