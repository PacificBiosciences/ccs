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

#include <ConsensusCore/Quiver/ReadScorer.hpp>

#include <ConsensusCore/Align/PairwiseAlignment.hpp>
#include <ConsensusCore/Quiver/QuiverConfig.hpp>
#include <ConsensusCore/Quiver/QvEvaluator.hpp>
#include <ConsensusCore/Quiver/SseRecursor.hpp>

#include <iostream>
#include <string>


using std::string;
using std::cout;
using std::endl;

namespace ConsensusCore
{
    ReadScorer::ReadScorer(QuiverConfig& config)
        : _quiverConfig(config)
    {}

    float ReadScorer::Score(const string& tpl, const QvRead& read) const
        throw(AlphaBetaMismatchException)
    {
        int I, J;
        SparseSseQvRecursor r(_quiverConfig.MovesAvailable, _quiverConfig.Banding);
        QvEvaluator e(read, tpl, _quiverConfig.QvParams);

        I = read.Length();
        J = tpl.length();
        SparseMatrixF alpha(I+1, J+1), beta(I+1, J+1);
        r.FillAlphaBeta(e, alpha, beta);

        return beta(0, 0);
    }

    const PairwiseAlignment*
    ReadScorer::Align(const string& tpl, const QvRead& read) const
        throw(AlphaBetaMismatchException)
    {
        int I, J;
        SparseSseQvRecursor r(_quiverConfig.MovesAvailable, _quiverConfig.Banding);
        QvEvaluator e(read, tpl, _quiverConfig.QvParams);

        I = read.Length();
        J = tpl.length();
        SparseMatrixF alpha(I+1, J+1), beta(I+1, J+1);
        r.FillAlphaBeta(e, alpha, beta);

        return r.Alignment(e, alpha);
    }

    const SparseMatrixF*
    ReadScorer::Alpha(const string& tpl, const QvRead& read) const
        throw(AlphaBetaMismatchException)
    {
        int I, J;
        SparseSseQvRecursor r(_quiverConfig.MovesAvailable, _quiverConfig.Banding);
        QvEvaluator e(read, tpl, _quiverConfig.QvParams);

        I = read.Length();
        J = tpl.length();
        SparseMatrixF *alpha = new SparseMatrixF(I+1, J+1);
        SparseMatrixF *beta  = new SparseMatrixF(I+1, J+1);
        r.FillAlphaBeta(e, *alpha, *beta);
        delete beta;
        return alpha;
    }

    const SparseMatrixF*
    ReadScorer::Beta(const string& tpl, const QvRead& read) const
        throw(AlphaBetaMismatchException)
    {
        int I, J;
        SparseSseQvRecursor r(_quiverConfig.MovesAvailable, _quiverConfig.Banding);
        QvEvaluator e(read, tpl, _quiverConfig.QvParams);

        I = read.Length();
        J = tpl.length();
        SparseMatrixF *alpha = new SparseMatrixF(I+1, J+1);
        SparseMatrixF *beta  = new SparseMatrixF(I+1, J+1);
        r.FillAlphaBeta(e, *alpha, *beta);
        delete alpha;
        return beta;
    }
}
