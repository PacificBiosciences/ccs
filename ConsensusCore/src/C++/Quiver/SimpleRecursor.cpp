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

#include <ConsensusCore/Quiver/SimpleRecursor.hpp>

#include <ConsensusCore/Edna/EdnaEvaluator.hpp>
#include <ConsensusCore/Matrix/DenseMatrix.hpp>
#include <ConsensusCore/Matrix/SparseMatrix.hpp>
#include <ConsensusCore/Quiver/detail/Combiner.hpp>
#include <ConsensusCore/Quiver/detail/RecursorBase.hpp>
#include <ConsensusCore/Quiver/QvEvaluator.hpp>
#include <ConsensusCore/Interval.hpp>
#include <ConsensusCore/Utils.hpp>

#include <algorithm>
#include <boost/tuple/tuple.hpp>
#include <climits>
#include <utility>

using std::min;
using std::max;

#define NEG_INF -FLT_MAX

namespace ConsensusCore {

    template<typename M, typename E, typename C>
    void
    SimpleRecursor<M, E, C>::FillAlpha(const E& e, const M& guide, M& alpha) const
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

            int i;
            float score = NEG_INF;
            float thresholdScore = NEG_INF;
            float maxScore = NEG_INF;

            alpha.StartEditingColumn(j, hintBeginRow, hintEndRow);

            int beginRow = hintBeginRow, endRow;
            for (i = beginRow;
                 i < I + 1 && (score >= thresholdScore || i < requiredEndRow);
                 ++i)
            {
                float thisMoveScore;
                score = NEG_INF;

                // Start:
                if (i == 0 && j == 0)
                {
                    score = 0.0f;
                }

                // Incorporation:
                if (i > 0 && j > 0)
                {
                    thisMoveScore = alpha(i - 1, j - 1) + e.Inc(i - 1, j - 1);
                    score = C::Combine(score, thisMoveScore);
                }

                // Extra:
                if (i > 0)
                {
                    thisMoveScore = alpha(i - 1, j) + e.Extra(i - 1, j);
                    score = C::Combine(score, thisMoveScore);
                }

                // Delete:
                if (j > 0)
                {
                    thisMoveScore = alpha(i, j - 1) + e.Del(i, j - 1);
                    score = C::Combine(score, thisMoveScore);
                }

                // Merge:
                if ((this->movesAvailable_ & MERGE) && j > 1 && i > 0)
                {
                    thisMoveScore = alpha(i - 1, j - 2) + e.Merge(i - 1, j - 2);
                    score = C::Combine(score, thisMoveScore);
                }

                //  Save score
                alpha.Set(i, j, score);

                if (score > maxScore)
                {
                    maxScore = score;
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
    SimpleRecursor<M, E, C>::FillBeta(const E& e, const M& guide, M& beta) const
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

            beta.StartEditingColumn(j, hintBeginRow, hintEndRow);

            int i;
            float score = NEG_INF;
            float thresholdScore = NEG_INF;
            float maxScore = NEG_INF;

            int beginRow, endRow = hintEndRow;
            for (i = endRow - 1;
                 i >= 0 && (score >= thresholdScore || i >= requiredBeginRow);
                 --i)
            {
                float thisMoveScore;
                score = NEG_INF;

                // Start:
                if (i == I && j == J)
                {
                    score = 0.0f;
                }

                // Incorporation:
                if (i < I && j < J)
                {
                    thisMoveScore = beta(i + 1, j + 1) + e.Inc(i, j);
                    score = C::Combine(score, thisMoveScore);
                }

                // Extra:
                if (i < I)
                {
                    thisMoveScore = beta(i + 1, j) + e.Extra(i, j);
                    score = C::Combine(score, thisMoveScore);
                }

                // Delete:
                if (j < J)
                {
                    thisMoveScore = beta(i, j + 1) + e.Del(i, j);
                    score = C::Combine(score, thisMoveScore);
                }

                // Merge:
                if ((this->movesAvailable_ & MERGE) && j < J - 1 && i < I)
                {
                    thisMoveScore = beta(i + 1, j + 2) + e.Merge(i, j);
                    score = C::Combine(score, thisMoveScore);
                }

                //  Save score
                beta.Set(i, j, score);

                if (score > maxScore)
                {
                    maxScore = score;
                    thresholdScore = maxScore - this->bandingOptions_.ScoreDiff;
                }
            }

            beginRow = i + 1;
            beta.FinishEditingColumn(j, beginRow, endRow);

            // Now, revise the hints to tell the caller where the mass of the
            // distribution really lived in this column.
            hintBeginRow = beginRow;
            for (i = endRow;
                 i > beginRow && beta(i - 1, j) < thresholdScore;
                 --i);
            hintEndRow = i;
        }
    }

    /// Calculate the recursion score by "stitching" together partial
    /// alpha and beta matrices.  alphaColumn, betaColumn, and
    /// absoluteColumn all refer to the same logical position in the
    /// template, but may have different values if, for instance,
    /// alpha here is a sub-range of the columns of the full alpha
    /// matrix.  Columns betaColumn and betaColumn + 1 of beta will be
    /// read; columns alphaColumn - 1 and alphaColumn - 2 of alpha
    /// will be read.
    template<typename M, typename E, typename C>
    float
    SimpleRecursor<M, E, C>::LinkAlphaBeta(const E& e,
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

        float v = NEG_INF, thisMoveScore;

        for (int i = usedBegin; i < usedEnd; i++)
        {
            if (i < I)
            {
                // Incorporate
                thisMoveScore = alpha(i, alphaColumn - 1) +
                                e.Inc(i, absoluteColumn - 1) +
                                beta(i + 1, betaColumn);
                v = C::Combine(v, thisMoveScore);

                // Merge (2 possible ways):
                thisMoveScore = alpha(i, alphaColumn - 2) +
                                e.Merge(i, absoluteColumn - 2) +
                                beta(i + 1, betaColumn);
                v = C::Combine(v, thisMoveScore);

                thisMoveScore = alpha(i, alphaColumn - 1) +
                                e.Merge(i, absoluteColumn - 1) +
                                beta(i + 1, betaColumn + 1);
                v = C::Combine(v, thisMoveScore);
            }

            // Delete:
            thisMoveScore = alpha(i, alphaColumn - 1) +
                            e.Del(i, absoluteColumn - 1) +
                            beta(i, betaColumn);
            v = C::Combine(v, thisMoveScore);
        }

        return v;
    }


    //
    // Reads: alpha(:, (beginColumn-2)..)
    //
    template<typename M, typename E, typename C>
    void
    SimpleRecursor<M, E, C>::ExtendAlpha(const E& e,
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
            float score;

            for (i = beginRow; i < endRow; i++)
            {
                float thisMoveScore;
                score = NEG_INF;

                // Incorporation:
                if (i > 0 && j > 0)
                {
                    float prev = extCol == 0 ?
                            alpha(i - 1, j - 1) :
                            ext(i - 1, extCol - 1);
                    thisMoveScore = prev + e.Inc(i - 1, j - 1);
                    score = C::Combine(score, thisMoveScore);
                }

                // Extra:
                if (i > 0)
                {
                    thisMoveScore = ext(i - 1, extCol) + e.Extra(i - 1, j);
                    score = C::Combine(score, thisMoveScore);
                }

                // Delete:
                if (j > 0)
                {
                    float prev = extCol == 0 ?
                            alpha(i, j - 1) :
                            ext(i, extCol - 1);
                    thisMoveScore = prev + e.Del(i, j - 1);
                    score = C::Combine(score, thisMoveScore);
                }

                // FIXME: is the merge code below incorrect for numExtColumns > 2?
                // Merge:
                if ((this->movesAvailable_ & MERGE) && j > 1 && i > 0)
                {
                    float prev = alpha(i - 1, j - 2);
                    thisMoveScore = prev + e.Merge(i - 1, j - 2);
                    score = C::Combine(score, thisMoveScore);
                }

                ext.Set(i, extCol, score);
            }
            assert (i == endRow);
            ext.FinishEditingColumn(extCol, beginRow, endRow);
        }
    }


    // Semantic: After ExtendBeta(B, j), we have
    //    ext(:, numExtColumns-1) = B'(:,j)
    //    ext(:, numExtColumns-2) = B'(:,j-1) ...
    //
    // Note: lastColumn is the numerically largest column number that
    // will be filled, but it is filled first since beta fill is done
    // backwards.
    //
    // Accesses B(:, ..(j+2))
    template<typename M, typename E, typename C>
    void
    SimpleRecursor<M, E, C>::ExtendBeta(const E& e,
                                        const M& beta, int lastColumn,
                                        M& ext, int numExtColumns,
                                        int lengthDiff) const
    {
        int I = beta.Rows() - 1;
        int J = beta.Columns() - 1;

        int lastExtColumn = numExtColumns - 1;

        assert(beta.Rows() == I + 1 &&
               ext.Rows() == I + 1);

        // The new template may not be the same length as the old template.
        // Just make sure that we have anough room to fill out the extend buffer
        assert(lastColumn + 2 <= J);
        assert(lastColumn >= 0);
        assert(ext.Columns() >= numExtColumns);

        for (int j = lastColumn; j > lastColumn - numExtColumns; j--)
        {
            int jp = j + lengthDiff;
            int extCol = lastExtColumn - (lastColumn - j);
            int beginRow, endRow;

            if (j < 0)
            {
                beginRow = 0;
                endRow = beta.UsedRowRange(0).End;
            }
            else
            {
                boost::tie(beginRow, endRow) = beta.UsedRowRange(j);
            }

            ext.StartEditingColumn(extCol, beginRow, endRow);

            int i;
            float score;

            for (i = endRow - 1;
                 i >= beginRow;
                 i--)
            {
                float thisMoveScore;
                score = NEG_INF;

                // Incorporation:
                if (i < I && j < J)
                {
                    float prev = (extCol == lastExtColumn) ?
                        beta(i + 1, j + 1) :
                        ext(i + 1, extCol + 1);
                    thisMoveScore = prev + e.Inc(i, jp);
                    score = C::Combine(score, thisMoveScore);
                }

                // Extra:
                if (i < I)
                {
                    thisMoveScore = ext(i + 1, extCol) + e.Extra(i, jp);
                    score = C::Combine(score, thisMoveScore);
                }

                // Delete:
                if (j < J)
                {
                    float prev = (extCol == lastExtColumn) ?
                        beta(i, j + 1) :
                        ext(i, extCol + 1);
                    thisMoveScore = prev + e.Del(i, jp);
                    score = C::Combine(score, thisMoveScore);
                }

                // FIXME: is the merge code below incorrect for numExtColumns > 2?
                // Merge:
                if ((this->movesAvailable_ & MERGE) && j < J - 1 && i < I)
                {
                    thisMoveScore = beta(i + 1, j + 2) + e.Merge(i, jp);
                    score = C::Combine(score, thisMoveScore);
                }

                ext.Set(i, extCol, score);
            }
            ext.FinishEditingColumn(extCol, beginRow, endRow);
        }
    }


    template<typename M, typename E, typename C>
    SimpleRecursor<M, E, C>::SimpleRecursor(int movesAvailable, const BandingOptions& banding)
        : detail::RecursorBase<M, E, C>(movesAvailable, banding)
    {}


    template class SimpleRecursor<DenseMatrixF,  QvEvaluator, detail::ViterbiCombiner>;
    template class SimpleRecursor<SparseMatrixF, QvEvaluator, detail::ViterbiCombiner>;
    template class SimpleRecursor<SparseMatrixF, QvEvaluator, detail::SumProductCombiner>;
    template class SimpleRecursor<SparseMatrixF, EdnaEvaluator, detail::SumProductCombiner>;
}
