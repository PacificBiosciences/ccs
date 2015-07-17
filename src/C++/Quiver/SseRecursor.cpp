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

#include <ConsensusCore/Quiver/SseRecursor.hpp>

#include <ConsensusCore/Interval.hpp>
#include <ConsensusCore/Utils.hpp>
#include <ConsensusCore/Edna/EdnaEvaluator.hpp>
#include <ConsensusCore/Matrix/DenseMatrix.hpp>
#include <ConsensusCore/Matrix/SparseMatrix.hpp>
#include <ConsensusCore/Quiver/detail/Combiner.hpp>
#include <ConsensusCore/Quiver/QvEvaluator.hpp>
#include <ConsensusCore/Quiver/SimpleRecursor.hpp>

#include <algorithm>
#include <boost/tuple/tuple.hpp>
#include <climits>
#include <numeric>
#include <utility>


using std::max;
using std::min;

#define NEG_INF   -FLT_MAX
#define POS_INF    FLT_MAX
#define NEG_INF_4  (Zero4<lfloat>())

namespace ConsensusCore {

// MSVC dones't appear to have this defined.
#ifdef _MSC_VER
    static inline __m128 operator+(const __m128 a, const __m128 b) {
        return _mm_add_ps(a, b);
    }
#endif


    template<typename M, typename E, typename C>
    void
    SseRecursor<M, E, C>::FillAlpha(const E& e, const M& guide, M& alpha) const
    {
        int I = e.ReadLength();
        int J = e.TemplateLength();

        assert(alpha.Rows() == I + 1 && alpha.Columns() == J + 1);
        assert(guide.IsNull() ||
               (guide.Rows() == alpha.Rows() && guide.Columns() == alpha.Columns()));

        int hintBeginRow = 0, hintEndRow = 0;

        for (int j = 0; j <= J; ++j)
        {
            this->RangeGuide(j, guide, alpha, &hintBeginRow, &hintEndRow);

            int requiredEndRow = min(I + 1, hintEndRow);

            float score = NEG_INF;
            float thresholdScore = NEG_INF;
            float maxScore = NEG_INF;

            alpha.StartEditingColumn(j, hintBeginRow, hintEndRow);

            int i;
            int beginRow = hintBeginRow, endRow;
            // Handle beginning rows non-SSE.  Must handle row 0 this
            // way (if row 0 is to be filled), and must terminate with
            // (I - i + 1) divisible by 4, so that the SSE loop can
            // run safely to the end.  Banding optimizations not applied
            // here.
            //
            // TODO(dalexander): we could also handle the first two columns this way,
            // and then we could remove all the conditionals from the SSE loop.
            // Profile first.
            //
            for (i = beginRow;
                 (i == 0 || (I - i + 1) % 4 != 0) && i <= I;
                 i++)
            {
                score = NEG_INF;

                // Start:
                if (i == 0 && j == 0)
                {
                    score = 0.0f;
                }
                // Inc
                if (i > 0 && j > 0)
                {
                    score = C::Combine(score, alpha(i - 1, j - 1) + e.Inc(i - 1, j - 1));
                }
                // Merge
                if ((this->movesAvailable_ & MERGE) && (i > 0 && j > 1))
                {
                    score = C::Combine(score, alpha(i - 1, j - 2) + e.Merge(i - 1, j - 2));
                }
                // Delete
                if (j > 0)
                {
                    score = C::Combine(score, alpha(i, j - 1) + e.Del(i, j - 1));
                }
                // Extra
                if (i > 0)
                {
                    score = C::Combine(score, alpha(i - 1, j) + e.Extra(i - 1, j));
                }
                alpha.Set(i, j, score);

                if (score > maxScore)
                {
                    maxScore = score;
                    thresholdScore = maxScore - this->bandingOptions_.ScoreDiff;
                }
            }
            //
            // Main SSE loop
            //
            assert(i > 0);
            for (;
                 i <= I && (score >= thresholdScore || i < requiredEndRow);
                 i += 4)
            {
                __m128 score4 = NEG_INF_4;
                // Incorporation:
                if (j > 0)
                {
                    score4 = C::Combine4(score4, alpha.Get4(i - 1, j - 1) + e.Inc4(i - 1, j - 1));
                }
                // Merge
                if ((this->movesAvailable_ & MERGE) && j >= 2)
                {
                    score4 = C::Combine4(score4, alpha.Get4(i - 1, j - 2) + e.Merge4(i - 1, j - 2));
                }
                // Deletion:
                if (j > 0)
                {
                    score4 = C::Combine4(score4, alpha.Get4(i, j - 1) + e.Del4(i, j - 1));
                }

                //
                // Extra (non-SSE cascade)
                //
                float insScores4_[4], score5_[5];

                __m128 insScores4 = e.Extra4(i - 1, j);
                _mm_storeu_ps(insScores4_, insScores4);

                score5_[0] = alpha.Get(i - 1, j);
                _mm_storeu_ps(&score5_[1], score4);

                for (int ii = 1; ii < 5; ii++)
                {
                    float v = C::Combine(score5_[ii], score5_[ii - 1] + insScores4_[ii - 1]);
                    score5_[ii] = v;
                }
                score4 = _mm_loadu_ps(&score5_[1]);
                alpha.Set4(i, j, score4);

                // Update score, potentialNewMax
                float potentialNewMax = *std::max_element(score5_ + 1, score5_ + 5);
                score = *std::min_element(score5_ + 1, score5_ + 5);

                if (potentialNewMax > maxScore)
                {
                    maxScore = potentialNewMax;
                    thresholdScore = maxScore - this->bandingOptions_.ScoreDiff;
                }
            }

            endRow = i;
            alpha.FinishEditingColumn(j, beginRow, endRow);

            // Now, revise the hints to tell the caller where the mass of the
            // distribution really lived in this column.
            hintEndRow = endRow;
            for (i = beginRow; i < endRow && alpha(i, j) < thresholdScore; ++i);
            hintBeginRow = i;
        }
    }


    template<typename M, typename E, typename C>
    void
    SseRecursor<M, E, C>::FillBeta(const E& e, const M& guide, M& beta) const
    {
        int I = e.ReadLength();
        int J = e.TemplateLength();

        assert(beta.Rows() == I + 1 && beta.Columns() == J + 1);
        assert(guide.IsNull() ||
               (guide.Rows() == beta.Rows() && guide.Columns() == beta.Columns()));

        int hintBeginRow = I + 1, hintEndRow = I + 1;

        for (int j = J; j >= 0; --j)
        {
            this->RangeGuide(j, guide, beta, &hintBeginRow, &hintEndRow);

            int requiredBeginRow = max(0, hintBeginRow);

            float score = NEG_INF;
            float thresholdScore = NEG_INF;
            float maxScore = NEG_INF;

            beta.StartEditingColumn(j, hintBeginRow, hintEndRow);
            //
            // See comment in FillAlpha---we are doing the same thing here.
            // An initial non-SSE loop, terminating when a multiple of 4
            // rows remain.
            //
            int i, beginRow, endRow = hintEndRow;
            for (i = endRow - 1;
                 (i == I || (i + 1) % 4 != 0) && i >= 0;
                 i--)
            {
                score = NEG_INF;

                // Start:
                if (i == I && j == J)
                {
                    score = 0.0f;
                }
                // Inc
                if (i < I && j < J)
                {
                    score = C::Combine(score, beta(i + 1, j + 1) + e.Inc(i, j));
                }
                // Merge
                if ((this->movesAvailable_ & MERGE) && j < J - 1 && i < I)
                {
                    score = C::Combine(score, beta(i + 1, j + 2) + e.Merge(i, j));
                }
                // Delete
                if (j < J)
                {
                    score = C::Combine(score, beta(i, j + 1) + e.Del(i, j));
                }
                // Extra
                if (i < I)
                {
                    score = C::Combine(score, beta(i + 1, j) + e.Extra(i, j));
                }

                beta.Set(i, j, score);

                if (score > maxScore)
                {
                    maxScore = score;
                    thresholdScore = maxScore - this->bandingOptions_.ScoreDiff;
                }
            }
            //
            // SSE loop
            //
            i = i - 3;
            for (;
                 i >= 0 && (score >= thresholdScore || i >= requiredBeginRow);
                 i -= 4)
            {
                __m128 score4 = NEG_INF_4;

                // Incorporation:
                if (i < I && j < J)
                {
                    score4 = C::Combine4(score4, beta.Get4(i + 1, j + 1) + e.Inc4(i, j));
                }
                // Merge
                if ((this->movesAvailable_ & MERGE) && j < J - 1 && i < I)
                {
                    score4 = C::Combine4(score4, beta.Get4(i + 1, j + 2) + e.Merge4(i, j));
                }
                // Deletion:
                if (j < J)
                {
                    score4 = C::Combine4(score4, beta.Get4(i, j + 1) + e.Del4(i, j));
                }

                //
                // Extra (non-SSE cascade)
                //
                float insScores4_[4], score5_[5];

                __m128 insScores4 = e.Extra4(i, j);
                _mm_storeu_ps(insScores4_, insScores4);

                score5_[4] = beta.Get(i + 4, j);
                _mm_storeu_ps(score5_, score4);

                for (int ii = 3; ii >= 0; ii--)
                {
                    float v = C::Combine(score5_[ii], score5_[ii + 1] + insScores4_[ii]);
                    score5_[ii] = v;
                }
                score4 = _mm_loadu_ps(score5_);
                beta.Set4(i, j, score4);

                // Update score, potentialNewMax
                float potentialNewMax = *std::max_element(score5_, score5_ + 4);
                score = *std::min_element(score5_, score5_ + 4);

                if (potentialNewMax > maxScore)
                {
                    maxScore = potentialNewMax;
                    thresholdScore = maxScore - this->bandingOptions_.ScoreDiff;
                }
            }

            beginRow = i + 4;
            beta.FinishEditingColumn(j, beginRow, endRow);

            // Now, revise the hints to tell the caller where the mass of the
            // distribution really lived in this column.
            hintBeginRow = beginRow;
            for (i = endRow;
                 i > beginRow && beta(i - 1, j) < thresholdScore;
                 i--);
            hintEndRow = i;
        }
    }

    template<typename M, typename E, typename C>
    INLINE_CALLEES float
    SseRecursor<M, E, C>::LinkAlphaBeta(const E& e,
                                        const M& alpha, int alphaColumn,
                                        const M& beta, int betaColumn,
                                        int absoluteColumn) const
    {
        const int I = e.ReadLength();

        assert(alphaColumn > 1 && absoluteColumn > 1);
        assert(absoluteColumn < e.TemplateLength());

        int usedBegin, usedEnd;
        boost::tie(usedBegin, usedEnd) = \
            RangeUnion(alpha.UsedRowRange(alphaColumn - 2),
                       alpha.UsedRowRange(alphaColumn - 1),
                       beta.UsedRowRange(betaColumn),
                       beta.UsedRowRange(betaColumn + 1));

        float v = NEG_INF;
        __m128 v4 = NEG_INF_4;

        // SSE loop
        int i;
        for (i = usedBegin; i < usedEnd - 4; i += 4)
        {
            // Incorporate
            v4 = C::Combine4(v4, alpha.Get4(i, alphaColumn - 1) +
                                 e.Inc4(i, absoluteColumn - 1) +
                                 beta.Get4(i + 1, betaColumn));
            // Merge (2 possible ways):
            if (this->movesAvailable_ & MERGE)
            {
                v4 = C::Combine4(v4, alpha.Get4(i, alphaColumn - 2) +
                                     e.Merge4(i, absoluteColumn - 2) +
                                     beta.Get4(i + 1, betaColumn));
                v4 = C::Combine4(v4, alpha.Get4(i, alphaColumn - 1) +
                                     e.Merge4(i, absoluteColumn - 1) +
                                     beta.Get4(i + 1, betaColumn + 1));
            }
            // Delete
            v4 = C::Combine4(v4, alpha.Get4(i, alphaColumn - 1) +
                                 e.Del4(i, absoluteColumn - 1) +
                                 beta.Get4(i, betaColumn));
        }
        // Handle the remaining rows non-SSE
        for (; i < usedEnd; i++)
        {
            if (i < I)
            {
                // Incorporate
                v = C::Combine(v, alpha(i, alphaColumn - 1) +
                                  e.Inc(i, absoluteColumn - 1) +
                                  beta(i + 1, betaColumn));
                // Merge (2 possible ways):
                if (this->movesAvailable_ & MERGE)
                {
                    v = C::Combine(v, alpha(i, alphaColumn - 2) +
                                      e.Merge(i, absoluteColumn - 2) +
                                      beta(i + 1, betaColumn));
                    v = C::Combine(v, alpha(i, alphaColumn - 1) +
                                      e.Merge(i, absoluteColumn - 1) +
                                      beta(i + 1, betaColumn + 1));
                }
            }
            // Delete:
            v = C::Combine(v, alpha(i, alphaColumn - 1) +
                              e.Del(i, absoluteColumn - 1) +
                              beta(i, betaColumn));
        }
        // Combine v4 and v
        float v_array[5];
        _mm_storeu_ps(v_array, v4);
        v_array[4] = v;
        v = std::accumulate(v_array, v_array + 5, NEG_INF, C::Combine);
        return v;
    }

    template<typename M, typename E, typename C>
    INLINE_CALLEES void
    SseRecursor<M, E, C>::ExtendAlpha(const E& e,
                                      const M& alpha, int beginColumn,
                                      M& ext, int numExtColumns) const
    {
        assert(numExtColumns >= 2);
        assert(alpha.Rows() == e.ReadLength() + 1 &&
               ext.Rows() == e.ReadLength() + 1);

        // The new template may not be the same length as the old template.
        // Just make sure that we have anough room to fill out the extend buffer
        assert(beginColumn + 1 < e.TemplateLength() + 1);
        assert(ext.Columns() >= numExtColumns);
        assert(beginColumn >= 2);

        for (int extCol = 0; extCol < numExtColumns; extCol++)
        {
            int j = beginColumn + extCol;
            int beginRow, endRow;

            //
            // If this extend is contained within the column bounds of
            // the original alpha, we use the row range that was
            // previously determined.  Otherwise start at alpha's last
            // UsedRow beginRow and go to the end.
            //
            if (j < alpha.Columns())
            {
                boost::tie(beginRow, endRow) = alpha.UsedRowRange(j);
            }
            else
            {
                beginRow = alpha.UsedRowRange(alpha.Columns() - 1).Begin;
                endRow = alpha.Rows();
            }

            ext.StartEditingColumn(extCol, beginRow, endRow);
            int i;
            // Handle the first rows non-SSE, leaving a multiple of 4
            // entries to be handed off to the SSE loop.  Need to always
            // handle at least row 0 this way, so that we don't have
            // to check for (i > 0) in the SSE loop.
            for (i = beginRow;
                 (i == 0 || (endRow - i) % 4 != 0) && i < endRow;
                 i++)
            {
                float prev, score = NEG_INF;
                if (i > 0)
                {
                    // Inc
                    prev = (extCol == 0 ?
                                alpha(i - 1, j - 1) :
                                ext(i - 1, extCol - 1));
                    score = C::Combine(score, prev + e.Inc(i - 1, j - 1));

                    // Extra
                    prev = ext(i - 1, extCol);
                    score = C::Combine(score, prev + e.Extra(i - 1, j));

                    // Merge
                    if (this->movesAvailable_ & MERGE)
                    {
                        prev = alpha(i - 1, j - 2);
                        score = C::Combine(score, prev + e.Merge(i - 1, j - 2));
                    }
                }
                // Delete
                prev = (extCol == 0 ?
                            alpha(i, j - 1) :
                            ext(i, extCol - 1));
                score = C::Combine(score, prev + e.Del(i, j - 1));
                ext.Set(i, extCol, score);
            }
            for (; i < endRow - 3; i += 4)
            {
                __m128 prev4, score4 = NEG_INF_4;

                // Incorporation:
                prev4 = (extCol == 0 ?
                            alpha.Get4(i - 1, j - 1) :
                            ext.Get4(i - 1, extCol - 1));
                score4 = C::Combine4(score4, prev4 + e.Inc4(i - 1, j - 1));

                // Merge
                if ((this->movesAvailable_ & MERGE) && j >= 2)
                {
                    prev4 = alpha.Get4(i - 1, j - 2);
                    score4 = C::Combine4(score4, prev4 + e.Merge4(i - 1, j - 2));
                }

                // Deletion:
                prev4 = (extCol == 0 ?
                            alpha.Get4(i, j - 1) :
                            ext.Get4(i, extCol - 1));
                score4 = C::Combine4(score4, prev4 + e.Del4(i, j - 1));

                // Extras:
                float insScores4_[4], score5_[5];

                __m128 insScores4 = e.Extra4(i - 1, j);
                _mm_storeu_ps(insScores4_, insScores4);

                score5_[0] = ext.Get(i - 1, extCol);
                _mm_storeu_ps(&score5_[1], score4);

                for (int ii = 1; ii < 5; ii++)
                {
                    float v = C::Combine(score5_[ii], score5_[ii - 1] + insScores4_[ii - 1]);
                    score5_[ii] = v;
                }
                score4 = _mm_loadu_ps(&score5_[1]);
                ext.Set4(i, extCol, score4);
            }
            assert (i == endRow);

            ext.FinishEditingColumn(extCol, beginRow, endRow);
        }
    }


    template<typename M, typename E, typename C>
    void
    SseRecursor<M, E, C>::ExtendBeta(const E& e,
                                     const M& beta, int endColumn,
                                     M& ext, int numExtColumns,
                                     int lengthDiff) const
    {
        simpleRecursor_.ExtendBeta(e, beta, endColumn, ext, numExtColumns, lengthDiff);
    }



    template<typename M, typename E, typename C>
    SseRecursor<M, E, C>::SseRecursor(int movesAvailable, const BandingOptions& banding)
        : detail::RecursorBase<M, E, C>(movesAvailable, banding),
          simpleRecursor_(movesAvailable, banding)
    {}

    template class SseRecursor<DenseMatrixF,  QvEvaluator, detail::ViterbiCombiner>;
    template class SseRecursor<SparseMatrixF, QvEvaluator, detail::ViterbiCombiner>;
    template class SseRecursor<SparseMatrixF, QvEvaluator, detail::SumProductCombiner>;
    template class SseRecursor<SparseMatrixF, EdnaEvaluator, detail::SumProductCombiner>;
}
