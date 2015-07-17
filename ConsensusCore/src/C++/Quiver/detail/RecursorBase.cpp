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

#include <ConsensusCore/Quiver/detail/RecursorBase.hpp>

#include <ConsensusCore/Align/PairwiseAlignment.hpp>
#include <ConsensusCore/LValue.hpp>
#include <ConsensusCore/Logging.hpp>
#include <ConsensusCore/Matrix/DenseMatrix.hpp>
#include <ConsensusCore/Matrix/SparseMatrix.hpp>
#include <ConsensusCore/Quiver/detail/Combiner.hpp>
#include <ConsensusCore/Edna/EdnaEvaluator.hpp>
#include <ConsensusCore/Quiver/QuiverConfig.hpp>
#include <ConsensusCore/Quiver/QvEvaluator.hpp>
#include <ConsensusCore/Types.hpp>
#include <ConsensusCore/Utils.hpp>

#include <algorithm>
#include <boost/type_traits.hpp>
#include <string>
#include <vector>


// TODO(dalexander): put these into a RecursorConfig struct
#define MAX_FLIP_FLOPS                  5
#define ALPHA_BETA_MISMATCH_TOLERANCE   0.2
#define REBANDING_THRESHOLD             0.04

using std::max;
using std::min;


namespace ConsensusCore {
namespace detail {

    template<typename M, typename E, typename C>
    int
    RecursorBase<M, E, C>::FillAlphaBeta(const E& e, M& a, M& b) const
        throw(AlphaBetaMismatchException)
    {
        FillAlpha(e, M::Null(), a);
        FillBeta(e, a, b);

        int I = e.ReadLength();
        int J = e.TemplateLength();
        int flipflops = 0;
        int maxSize = static_cast<int>(0.5 + REBANDING_THRESHOLD * (I + 1) * (J + 1));

        // if we use too much space, do at least one more round
        // to take advantage of rebanding
        if (a.UsedEntries() >= maxSize ||
            b.UsedEntries() >= maxSize)
        {
            FillAlpha(e, b, a);
            FillBeta(e, a, b);
            FillAlpha(e, b, a);
            flipflops += 3;
        }

        while (fabs(a(I, J) - b(0, 0)) > ALPHA_BETA_MISMATCH_TOLERANCE
               && flipflops <= MAX_FLIP_FLOPS)
        {
            if (flipflops % 2 == 0)
            {
                FillAlpha(e, b, a);
            }
            else
            {
                FillBeta(e, a, b);
            }
            flipflops++;
        }

        if (fabs(a(I, J) - b(0, 0)) > ALPHA_BETA_MISMATCH_TOLERANCE)
        {
            LDEBUG << "Could not mate alpha, beta.  Read: "
                   << e.ReadName() << " Tpl: " << e.Template();
            throw AlphaBetaMismatchException();
        }

        return flipflops;
    }

    struct MoveSpec {
        Move MoveType;
        int ReadDelta;
        int ReferenceDelta;
    };

    template<typename M, typename E, typename C>
    const PairwiseAlignment*
    RecursorBase<M, E, C>::Alignment(const E& e, const M& a) const
    {
        if (!boost::is_same<C, ViterbiCombiner>::value)
        {
            ShouldNotReachHere();
        }

        int I, J;
        I = e.ReadLength();
        J = e.TemplateLength();

        // Matrix must be filled in before requesting traceback
        assert(a(I, J) != lfloat());

        int i = I;
        int j = J;
        float pathScore = 0;

        MoveSpec incMove   = { INCORPORATE, 1, 1 };
        MoveSpec delMove   = { DELETE,      0, 1 };
        MoveSpec extraMove = { EXTRA,       1, 0 };
        MoveSpec mergeMove = { MERGE,       1, 2 };
        std::vector<MoveSpec> moves;

        while (i > 0 || j > 0)
        {
            MoveSpec bestMove = { INVALID_MOVE, 0, 0 };
            float prevScore, moveScore, totalScore;
            float bestScore = lfloat();
            float bestMoveScore = lfloat();

            if (i > 0 && j > 0)
            {
                prevScore = a(i - 1, j - 1);
                moveScore = e.Inc(i - 1, j - 1);
                totalScore = prevScore + moveScore;
                if (totalScore > bestScore)
                {
                    // Incorporate (match or mismatch)
                    bestMove = incMove;
                    bestScore = totalScore;
                    bestMoveScore = moveScore;
                }
            }

            if (j > 0)
            {
                // Delete
                prevScore = a(i, j - 1);
                bool freeDelete = (!e.PinEnd() && i == I) || (!e.PinStart() && i == 0);
                moveScore = freeDelete ? 0 : e.Del(i, j - 1);
                totalScore = prevScore + moveScore;
                if (totalScore > bestScore)
                {
                    bestMove = delMove;
                    bestScore = totalScore;
                    bestMoveScore = moveScore;
                }
            }

            if (i > 0)
            {
                // Extra
                prevScore = a(i - 1, j);
                moveScore = e.Extra(i - 1, j);
                totalScore = prevScore + moveScore;
                if (totalScore > bestScore)
                {
                    bestMove = extraMove;
                    bestScore = totalScore;
                    bestMoveScore = moveScore;
                }
            }

            if ((movesAvailable_ & MERGE) && i > 0 && j > 1)
            {
                // Merge
                prevScore = a(i - 1, j - 2);
                moveScore = e.Merge(i - 1, j - 2);
                totalScore = prevScore + moveScore;
                if (totalScore > bestScore)
                {
                    bestMove = mergeMove;
                    bestScore = totalScore;
                    bestMoveScore = moveScore;
                }
            }
            assert(AlmostEqual(a(i, j), bestScore));
            assert(bestMove.MoveType != INVALID_MOVE);
            assert(bestMoveScore != lfloat());

            moves.push_back(bestMove);
            i -= bestMove.ReadDelta;
            j -= bestMove.ReferenceDelta;
            pathScore += bestMoveScore;
        }
        assert(i == 0 && j == 0);

        // Reverse moves
        std::reverse(moves.begin(), moves.end());

        // Replay moves and stringify.
        std::string target;
        std::string query;
        i = 0;
        j = 0;
        foreach (const MoveSpec& move, moves)
        {
            switch (move.MoveType) {
            case INCORPORATE:
                target     += e.Template()[j];
                query      += e.Basecalls()[i];
                break;
            case EXTRA:
                target     += '-';
                query      += e.Basecalls()[i];
                break;
            case DELETE:
                target     += e.Template()[j];
                query      += '-';
                break;
            case MERGE:
                target     += e.Template()[j];
                target     += e.Template()[j+1];
                query      += '-';
                query      += e.Basecalls()[i];
                break;
            case INVALID_MOVE:
                ShouldNotReachHere();
                break;
            default:
                ShouldNotReachHere();
            }
            i += move.ReadDelta;
            j += move.ReferenceDelta;
            assert(target.length() == query.length());
        }

        return new PairwiseAlignment(target, query);
    }

    template<typename M, typename E, typename C>
    RecursorBase<M, E, C>::RecursorBase(int movesAvailable, const BandingOptions& bandingOptions)
        : movesAvailable_(movesAvailable),
          bandingOptions_(bandingOptions)
    {}

    template<typename M, typename E, typename C>
    RecursorBase<M, E, C>::~RecursorBase()
    {}

    // template instantiation
    template class RecursorBase<DenseMatrixF, QvEvaluator, ViterbiCombiner>;
    template class RecursorBase<SparseMatrixF, QvEvaluator, ViterbiCombiner>;
    template class RecursorBase<SparseMatrixF, QvEvaluator, SumProductCombiner>;
    template class RecursorBase<SparseMatrixF, EdnaEvaluator, SumProductCombiner>;
}}
