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

// The POA "alignment matrix" is a set of alignment columns
// corresponding to each vertex in the graph that the read was aligned
// against.

#pragma once

#include <pacbio/consensus/poa/PoaGraph.h>

#include <boost/unordered_map.hpp>
#include <boost/utility.hpp>
#include <cfloat>

#include "BoostGraph.h"
#include "VectorL.h"

using boost::noncopyable;
using boost::unordered_map;

namespace PacBio {
namespace Consensus {
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

typedef unordered_map<VD, const AlignmentColumn*> AlignmentColumnMap;

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
    AlignMode mode_;
    float score_;
};
}
}
}  // PacBio::Consensus::detail
