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

// Author: David Alexander, Lance Hepler

#pragma once

#include <algorithm>
#include <climits>
#include <memory>
#include <utility>

#include <pacbio/consensus/Exceptions.h>
#include <pacbio/consensus/Read.h>
#include <pacbio/consensus/Template.h>

#include "matrix/ScaledMatrix.h"

namespace PacBio {
namespace Consensus {

// AbstractRecursor is in Template.h

// TODO(lhepler) comment about use of CRTP
template <typename Derived>
class Recursor : public AbstractRecursor
{
public:
    // \brief Construct a Recursor from a Template and a MappedRead,
    // The scoreDiff here is passed in negative logScale and converted
    // to the appropriate divisor.
    Recursor(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr,
             double scoreDiff = 12.5);

    /// \brief Fill the alpha and beta matrices.
    /// This routine will fill the alpha and beta matrices, ensuring
    /// that the score computed from the alpha and beta recursions are
    /// identical, refilling back-and-forth if necessary.
    size_t FillAlphaBeta(M& alpha, M& beta) const throw(AlphaBetaMismatch);

    /**
     Fill in the alpha matrix.  This matrix has the read run along the rows, and
     the template run along the columns.  The first row and column do not correspond
     to a template position.  Therefore the match represented at position (i,j)
     corresponds to a match between template positions (i+1, j+1).

     The alpha matrix is the "Forward" matrix used in the forward/backward
     algorithm.
     The i,j position of the matrix represents the probability of all paths up
     to the point where the ith read position and jth template have been
     "emitted."
     The matrix is calculated recursively by examining all possible transitions
     into (i,j), and calculating the probability we were in the previous state,
     times the probability of a transition into (i,j) times the probability of
     emitting the observation that corresponds to (i,j). All probabilities are
     calculated and stored as LOG values.

     Note that in doing this calculation, in order to work with di-nucleotide
     contexts, we require that the first and last transition be a match.  In other words the
     start and end of the read and template are "pinned" to each other.

     //TODO: Verify memory is initialized to 0!

     @param guide An object that helps inform how to select the size of "bands"
     for the banded algorithm used.  This is typically the beta matrix if we are
     "repopulating" the matrix.
     @param alpha The matrix to be filled.
     */

    void FillAlpha(const M& guide, M& alpha) const;

    /**
     Fill the Beta matrix, the backwards half of the forward-backward algorithm.
     This represents the probability that starting from the (i,j) state, the
     combined probability of transitioning out and following all paths through to the
     end. That is, we need to calculate transition from state and emit from next
     state for each

     In combination with the Alpha matrix, this allows us to calculate all paths
     that pass through the (i,j) element, as exp(Alpha(i,j) + Beta(i,j))

     All probabilities stored in the matrix are stored as NON-LOGGED
     probabilities.

     @param e The evaluator, such as QvEvaluator
     @param M the guide matrix for banding (this needs more documentation)
     @param beta The Beta matrix, stored as either a DenseMatrix or a
     SparseMatrix.
     */

    void FillBeta(const M& guide, M& beta) const;

    /// \brief Calculate the recursion score by "linking" partial alpha and/or
    ///        beta matrices.
    double LinkAlphaBeta(const M& alpha, size_t alphaColumn, const M& beta, size_t betaColumn,
                         size_t absoluteColumn) const;

    void ExtendAlpha(const M& alpha, size_t beginColumn, M& ext, size_t numExtColumns = 2) const;

    void ExtendBeta(const M& beta, size_t endColumn, M& ext, int lengthDiff = 0) const;

private:
    std::pair<size_t, size_t> RowRange(size_t j, const M& matrix) const;

    /// \brief Reband alpha and beta matrices.
    /// This routine will reband alpha and beta to the convex hull
    /// of the maximum path through each and the inputs for column j.
    bool RangeGuide(size_t j, const M& guide, const M& matrix, size_t* beginRow,
                    size_t* endRow) const;
    // The RangeGuide function determines the minimum score by dividing out scoreDiff_

private:
    std::vector<uint8_t> emissions_;
};

namespace {  // anonymous

typedef std::pair<size_t, size_t> Interval;

// TODO(dalexander): put these into a RecursorConfig struct
// TODO(anybody): Hmmm... not sure what the heck to do about these...
constexpr int MAX_FLIP_FLOPS = 5;
constexpr double ALPHA_BETA_MISMATCH_TOLERANCE = 0.001;
constexpr double REBANDING_THRESHOLD = 0.04;

constexpr uint8_t kDefaultBase = 0;  // corresponding to A, usually
constexpr TemplatePosition kDefaultTplPos = TemplatePosition{'A', kDefaultBase, 1, 0, 0, 0};

inline Interval RangeUnion(const Interval& range1, const Interval& range2)
{
    return Interval(std::min(range1.first, range2.first), std::max(range1.second, range2.second));
}

inline Interval RangeUnion(const Interval& range1, const Interval& range2, const Interval& range3,
                           const Interval& range4)
{
    return RangeUnion(RangeUnion(range1, range2), RangeUnion(range3, range4));
}

inline double Combine(const double a, const double b) { return a + b; }
}  // namespace anonymous

template <typename Derived>
void Recursor<Derived>::FillAlpha(const M& guide, M& alpha) const
{
    // We are pinning, so should never go all the way to the end of the
    // read/template
    // But our matrix indexing is one off the model/outcome indexing
    // so the match in (1,1) corresponds to a pairing between
    // Model[0]/Outcome[0]
    size_t I = read_.Length();
    size_t J = tpl_->Length();

    assert(alpha.Rows() == I + 1 && alpha.Columns() == J + 1);
    assert(guide.IsNull() || (guide.Rows() == alpha.Rows() && guide.Columns() == alpha.Columns()));

    // Initial condition, we always start with a match
    alpha.StartEditingColumn(0, 0, 1);
    alpha.Set(0, 0, 1.0);
    alpha.FinishEditingColumn(0, 0, 1);
    // End initial conditions

    size_t hintBeginRow = 1, hintEndRow = 1;
    auto prevTransProbs = kDefaultTplPos;
    uint8_t prevTplBase = prevTransProbs.Idx;

    for (int j = 1; j < J; ++j)  // Note due to offset with reads and otherwise, this is ugly-ish
    {
        // Load up the transition parameters for this context

        auto currTransProbs = (*tpl_)[j - 1];
        auto currTplBase = currTransProbs.Idx;
        this->RangeGuide(j, guide, alpha, &hintBeginRow, &hintEndRow);

        size_t i;
        double thresholdScore = 0.0;
        double maxScore = 0.0;
        double score = 0.0;
        alpha.StartEditingColumn(j, hintBeginRow, hintEndRow);

        auto nextTplBase = (*tpl_)[j].Idx;

        size_t beginRow = hintBeginRow, endRow;
        // Recursively calculate [Probability in last state] * [Probability
        // transition to new state] * [Probability of emission]
        for (i = beginRow; i < I && (score >= thresholdScore || i < hintEndRow); ++i) {
            // TODO: Terrible hack right now to emit this guy as teh IQV
            const uint8_t curReadEm = emissions_[i - 1];
            double thisMoveScore = 0.0;
            score = 0.0;
            // Match:
            /* Important!  Note that because we require the initial state to be
              a match,
               when i = 1 and j = 1 the match transition probability must be 1,
              since no other options
               are allowed.  Similarly, the probability for the match
              probability to the end base should be 1.

               Note that for the first "match" between a read and template, we
              have no choice but to
               hard code it to 1, as there is no defined transition probability
              for a dinucleotide context.

              ***********  EDGE_CONDITION ************
             */
            if (i > 0 && j > 0) {
                thisMoveScore =
                    alpha(i - 1, j - 1) * prevTransProbs.Match *
                    Derived::EmissionPr(MoveType::MATCH, curReadEm, prevTplBase, currTplBase);
                score = Combine(score, thisMoveScore);
            }

            if (i > 1) {
                // Branch, due to pinning, can't "insert" first or last read base
                thisMoveScore =
                    alpha(i - 1, j) * currTransProbs.Branch *
                    Derived::EmissionPr(MoveType::BRANCH, curReadEm, currTplBase, nextTplBase);
                score = Combine(score, thisMoveScore);

                // Stick
                thisMoveScore =
                    alpha(i - 1, j) * currTransProbs.Stick *
                    Derived::EmissionPr(MoveType::STICK, curReadEm, currTplBase, nextTplBase);
                score = Combine(score, thisMoveScore);
            }

            // Deletion, due to pinning, can't "delete" first or last template bp
            if (j > 1) {
                thisMoveScore = alpha(i, j - 1) * prevTransProbs.Deletion;
                score = Combine(score, thisMoveScore);
            }

            //  Save score
            alpha.Set(i, j, score);

            if (score > maxScore) {
                maxScore = score;
                thresholdScore = maxScore / scoreDiff_;
            }
        }
        endRow = i;
        prevTransProbs = currTransProbs;
        prevTplBase = currTplBase;
        // Now, revise the hints to tell the caller where the mass of the
        // distribution really lived in this column.
        hintEndRow = endRow;
        for (i = beginRow; i < endRow && alpha(i, j) < thresholdScore; ++i)
            ;
        hintBeginRow = i;

        // Don't rescale until we finish updating the hint.
        alpha.FinishEditingColumn(j, beginRow, endRow);
    }

    /* Now fill out the probability in the last pinned position.
     * We require that we end in a match.
     * search for the term EDGE_CONDITION to find a comment with more
     * information */
    {
        auto currTplBase = (*tpl_)[J - 1].Idx;
        assert(J < 2 || prevTplBase == (*tpl_)[J - 2].Idx);
        // end in the homopolymer state for now.
        auto likelihood =
            alpha(I - 1, J - 1) *
            Derived::EmissionPr(MoveType::MATCH, emissions_[I - 1], prevTplBase, currTplBase);
        alpha.StartEditingColumn(J, I, I + 1);
        alpha.Set(I, J, likelihood);
        alpha.FinishEditingColumn(J, I, I + 1);
    }
}

template <typename Derived>
void Recursor<Derived>::FillBeta(const M& guide, M& beta) const
{
    size_t I = read_.Length();
    size_t J = tpl_->Length();

    assert(beta.Rows() == I + 1 && beta.Columns() == J + 1);
    assert(guide.IsNull() || (guide.Rows() == beta.Rows() && guide.Columns() == beta.Columns()));

    // Setup initial condition, at the end we are one
    beta.StartEditingColumn(J, I, I + 1);
    beta.Set(I, J, 1.0);
    beta.FinishEditingColumn(J, I, I + 1);

    // Totally arbitray decision here...
    size_t hintBeginRow = I, hintEndRow = I;

    // Recursively calculate [Probability transition to next state] *
    // [Probability of emission at that state] * [Probability from that state]
    for (int j = J - 1; j > 0; --j) {
        const auto nextTplPos = (*tpl_)[j];
        const auto nextTplBase = nextTplPos.Idx;
        const auto currTransProbs = (*tpl_)[j - 1];

        this->RangeGuide(j, guide, beta, &hintBeginRow, &hintEndRow);

        beta.StartEditingColumn(j, hintBeginRow, hintEndRow);

        int i;
        double score = 0.0;
        double thresholdScore = 0.0;
        double maxScore = 0.0;

        int beginRow, endRow = hintEndRow;
        for (i = endRow - 1; i > 0 && (score >= thresholdScore || i >= hintBeginRow); --i) {
            const uint8_t nextReadEm = emissions_[i];
            double thisMoveScore = 0.0;
            score = 0.0;

            // Match
            if (i + 1 < I) {
                thisMoveScore = beta(i + 1, j + 1) * currTransProbs.Match *
                                Derived::EmissionPr(MoveType::MATCH, nextReadEm, currTransProbs.Idx,
                                                    nextTplBase);
                score = Combine(score, thisMoveScore);
            } else if (i + 1 == I && j + 1 == J) {
                thisMoveScore =
                    beta(i + 1, j + 1) * Derived::EmissionPr(MoveType::MATCH, nextReadEm,
                                                             currTransProbs.Idx, nextTplBase);
                score = Combine(score, thisMoveScore);
            }

            // Branch, can only transition to an insertion for the 2nd to last read base and before
            if (0 < i && i < I) {
                thisMoveScore = beta(i + 1, j) * currTransProbs.Branch *
                                Derived::EmissionPr(MoveType::BRANCH, nextReadEm,
                                                    currTransProbs.Idx, nextTplBase);
                score = Combine(score, thisMoveScore);

                // Stick, can only transition to an insertion for the 2nd to last read base
                thisMoveScore = beta(i + 1, j) * currTransProbs.Stick *
                                Derived::EmissionPr(MoveType::STICK, nextReadEm, currTransProbs.Idx,
                                                    nextTplBase);
                score = Combine(score, thisMoveScore);
            }

            // Deletion
            if (0 < j && j < J) {
                thisMoveScore = beta(i, j + 1) * currTransProbs.Deletion;
                score = Combine(score, thisMoveScore);
            }

            // Save score
            beta.Set(i, j, score);

            if (score > maxScore) {
                maxScore = score;
                thresholdScore = maxScore / scoreDiff_;
            }
        }

        beginRow = i + 1;
        // DumpBetaMatrix(beta);
        // Now, revise the hints to tell the caller where the mass of the
        // distribution really lived in this column.

        hintBeginRow = beginRow;
        for (i = endRow; i > beginRow && beta(i - 1, j) < thresholdScore; --i)
            ;
        hintEndRow = i;

        // Don't rescale until we update the hints
        beta.FinishEditingColumn(j, beginRow, endRow);
    }

    /* Now to fill the top row which must be a match
     * search for the term EDGE_CONDITION to find a comment with more
     * information */
    {
        beta.StartEditingColumn(0, 0, 1);
        auto match_emission_prob =
            Derived::EmissionPr(MoveType::MATCH, emissions_[0], kDefaultBase, (*tpl_)[0].Idx);
        beta.Set(0, 0, match_emission_prob * beta(1, 1));
        beta.FinishEditingColumn(0, 0, 1);
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
template <typename Derived>
double Recursor<Derived>::LinkAlphaBeta(const M& alpha, size_t alphaColumn, const M& beta,
                                        size_t betaColumn, size_t absoluteColumn) const
{
    const size_t I = read_.Length();

    assert(alphaColumn > 1 && absoluteColumn > 1);
    assert(absoluteColumn <= tpl_->Length());

    int usedBegin, usedEnd;
    std::tie(usedBegin, usedEnd) =
        RangeUnion(alpha.UsedRowRange(alphaColumn - 2), alpha.UsedRowRange(alphaColumn - 1),
                   beta.UsedRowRange(betaColumn), beta.UsedRowRange(betaColumn + 1));

    double v = 0.0, thisMoveScore;

    const auto currTplParams = (*tpl_)[absoluteColumn - 1];
    const auto prevTplParams = (*tpl_)[absoluteColumn - 2];

    for (int i = usedBegin; i < usedEnd; i++) {
        if (i < I) {
            const uint8_t readEm = emissions_[i];
            // Match
            thisMoveScore =
                alpha(i, alphaColumn - 1) * prevTplParams.Match *
                Derived::EmissionPr(MoveType::MATCH, readEm, prevTplParams.Idx, currTplParams.Idx) *
                beta(i + 1, betaColumn);
            v = Combine(v, thisMoveScore);
        }

        // Delete
        thisMoveScore = alpha(i, alphaColumn - 1) * prevTplParams.Deletion * beta(i, betaColumn);
        v = Combine(v, thisMoveScore);
    }

    return (std::log(v) + alpha.GetLogProdScales(0, alphaColumn) +
            beta.GetLogProdScales(betaColumn, beta.Columns()));
}

/**
 This method extends that Alpha matrix into a temporary matrix given by
 ext.  It extends the region [beginColumn, beginColumn + numExtColumns)

 Note that this method is used EXCLUSIVELY for testing mutations, and so
 we don't get the actual parameters and positions from the template, but we
 get them after a "virtual" mutation has been applied.

 All new data is placed in the extension matrix.  The guesses for start/end
 rows in the banding are determined by evaluating neighbors of each position.
 @param <#parameter#>
 @returns <#retval#>
 */
template <typename Derived>
void Recursor<Derived>::ExtendAlpha(const M& alpha, size_t beginColumn, M& ext,
                                    size_t numExtColumns) const
{
    assert(numExtColumns >= 2);  // We have to fill at least one
    assert(alpha.Rows() == read_.Length() + 1 &&
           ext.Rows() == read_.Length() + 1);  // The read never mutates

    // The new template may not be the same length as the old template.
    // Just make sure that we have anough room to fill out the extend buffer
    assert(beginColumn + 1 < tpl_->Length() + 1);
    assert(ext.Columns() >= numExtColumns);
    assert(beginColumn >= 2);
    // Due to pinning at the end, moves are only possible if less than these
    // positions.
    size_t maxLeftMovePossible = tpl_->Length();
    size_t maxDownMovePossible = read_.Length();

    // completely fill the rectangle bounded by the min and max
    size_t beginRow, endRow;
    std::tie(beginRow, endRow) = alpha.UsedRowRange(beginColumn);
    for (size_t j = 1; j + beginColumn < alpha.Columns() && j <= numExtColumns; ++j)
        endRow = std::max(alpha.UsedRowRange(j + beginColumn).second, endRow);

    for (size_t extCol = 0; extCol < numExtColumns; extCol++) {
        size_t j = beginColumn + extCol;

        ext.StartEditingColumn(extCol, beginRow, endRow);

        int i;
        double score = 0.0;

        // Grab values that will be useful for the whole column
        auto currTplParams = (*tpl_)[j - 1];
        auto currTplBase = currTplParams.Idx;
        TemplatePosition prevTplParams = kDefaultTplPos;
        if (j > 1) {
            prevTplParams = (*tpl_)[j - 2];
        }
        uint8_t nextTplBase;
        if (j != maxLeftMovePossible) {
            nextTplBase = (*tpl_)[j].Idx;
        }

        for (i = beginRow; i < endRow; i++) {
            const uint8_t currReadEm = emissions_[i - 1];
            double thisMoveScore = 0.0;

            // Match
            if (i > 0 && j > 0) {
                double prev = extCol == 0 ? alpha(i - 1, j - 1) : ext(i - 1, extCol - 1);
                if (i < maxDownMovePossible && j < maxLeftMovePossible) {
                    thisMoveScore = prev * prevTplParams.Match *
                                    Derived::EmissionPr(MoveType::MATCH, currReadEm,
                                                        prevTplParams.Idx, currTplParams.Idx);
                } else if (i == maxDownMovePossible && j == maxLeftMovePossible) {
                    thisMoveScore =
                        prev * Derived::EmissionPr(MoveType::MATCH, currReadEm, prevTplParams.Idx,
                                                   currTplParams.Idx);
                }
                score = thisMoveScore;
            }

            // Branch
            if (i > 1 && i < maxDownMovePossible && j != maxLeftMovePossible) {
                thisMoveScore =
                    ext(i - 1, extCol) * currTplParams.Branch *
                    Derived::EmissionPr(MoveType::BRANCH, currReadEm, currTplBase, nextTplBase);
                score = Combine(score, thisMoveScore);
            }

            // Stick
            if (i > 1 && i < maxDownMovePossible && j != maxLeftMovePossible) {
                thisMoveScore =
                    ext(i - 1, extCol) * currTplParams.Stick *
                    Derived::EmissionPr(MoveType::STICK, currReadEm, currTplBase, nextTplBase);
                score = Combine(score, thisMoveScore);
            }

            // Delete
            if (j > 1 && j < maxLeftMovePossible && i != maxDownMovePossible) {
                double prev = extCol == 0 ? alpha(i, j - 1) : ext(i, extCol - 1);
                thisMoveScore = prev * prevTplParams.Deletion;
                score = Combine(score, thisMoveScore);
            }
            ext.Set(i, extCol, score);
        }
        assert(i == endRow);
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

/* Note this is a very confusing routine in order to avoid recomputing and
   additional memory allocations.  This routine tries to stick on a beta matrix
   to the original and back trace to the 0,0 position of this extension matrix.
   matrix from the original.  Note that the original beta
   matrix is indexed by the original template positions, while the template
   bases and parameters are now indexed according the the "virtual" template
   to which mutations have been applied.
 */
// @param lastColumn - Where we
template <typename Derived>
void Recursor<Derived>::ExtendBeta(const M& beta, size_t lastColumn, M& ext, int lengthDiff) const
{
    size_t I = read_.Length();
    size_t J = tpl_->Length();

    // How far back do we have to go until we are at the zero (first) column?
    // we always go all the way back.
    int numExtColumns = 1 + lengthDiff + lastColumn;
    int firstColumn = 0 - lengthDiff;
    int lastExtColumn = numExtColumns - 1;

    // The new template may not be the same length as the old template.
    // Just make sure that we have enough room to fill out the extend buffer
    assert(lastColumn + 1 <= J);
    assert(lastColumn < 4);  // Since we are only testing mutations of size 1,
                             // and the check prior for a beginning mutation is
                             // < 3, max = 2 + 1 = 3
    assert(lastColumn >= 0);
    assert(ext.Columns() >= numExtColumns);
    assert(beta.Rows() == I + 1 && ext.Rows() == I + 1);
    assert(abs(lengthDiff) < 2);

    // completely fill the rectangle bounded by the min and max
    int beginRow, endRow;
    if (lastColumn + 1 < beta.Columns())
        std::tie(beginRow, endRow) = beta.UsedRowRange(lastColumn + 1);
    else
        std::tie(beginRow, endRow) = beta.UsedRowRange(lastColumn);
    for (int j = 0; j <= lastColumn && j <= numExtColumns; ++j)
        beginRow = std::min(static_cast<int>(beta.UsedRowRange(lastColumn - j).first), beginRow);

    for (int j = lastColumn; j + numExtColumns > lastColumn; j--) {
        /* Convert from old template to new template coordinates.
           lengthDiff will be 0 for substitution, -1 for deletion and +1 for
           insertion
         */
        int jp = j + lengthDiff;
        // What is the current extension column we are adding data into.
        int extCol = lastExtColumn - (lastColumn - j);

        ext.StartEditingColumn(extCol, beginRow, endRow);

        // Load up useful values referenced throughout the column.
        auto nextTplParams = (*tpl_)[jp];
        auto nextTplBase = nextTplParams.Idx;

        TemplatePosition currTplParams = kDefaultTplPos;
        if (jp > 0) currTplParams = (*tpl_)[jp - 1];

        for (int i = endRow - 1; i >= beginRow; i--) {
            uint8_t nextReadEm = 4;  // 'N'
            if (i < I) {
                nextReadEm = emissions_[i];
            }
            double thisMoveScore = 0.0;
            double score = 0.0;

            // Match
            // TODO: Remove these checks, we should always be on the left side
            // of the matrix....
            if (i < I && j < J) {
                double next =
                    (extCol == lastExtColumn) ? beta(i + 1, j + 1) : ext(i + 1, extCol + 1);

                // First and last have to start with an emission
                // TODO: So ugly we need to clean this up!
                // All these checks should be reorganized, redundand subexpressions
                // combined.
                if (j > firstColumn && i > 0) {
                    thisMoveScore = next * currTplParams.Match *
                                    Derived::EmissionPr(MoveType::MATCH, nextReadEm,
                                                        currTplParams.Idx, nextTplBase);
                    score = Combine(score, thisMoveScore);
                }
            }

            // Branch
            if (0 < i && i < I && firstColumn < j) {
                thisMoveScore = ext(i + 1, extCol) * currTplParams.Branch *
                                Derived::EmissionPr(MoveType::BRANCH, nextReadEm, currTplParams.Idx,
                                                    nextTplBase);
                score = Combine(score, thisMoveScore);
            }

            // Stick
            if (0 < i && i < I && firstColumn < j) {
                thisMoveScore = ext(i + 1, extCol) * currTplParams.Stick *
                                Derived::EmissionPr(MoveType::STICK, nextReadEm, currTplParams.Idx,
                                                    nextTplBase);
                score = Combine(score, thisMoveScore);
            }

            // Deletion
            if (0 < i && firstColumn < j && j < J) {
                double next = (extCol == lastExtColumn) ? beta(i, j + 1) : ext(i, extCol + 1);
                thisMoveScore = next * currTplParams.Deletion;
                score = Combine(score, thisMoveScore);
            }

            ext.Set(i, extCol, score);
        }
        ext.FinishEditingColumn(extCol, beginRow, endRow);
    }

    // fill out the (0, 0) entry of the matrix
    {
        ext.StartEditingColumn(0, 0, 1);
        const double match_trans_prob = (lastExtColumn == 0) ? beta(1, lastColumn + 1) : ext(1, 1);
        const double match_emission_prob =
            Derived::EmissionPr(MoveType::MATCH, emissions_[0], kDefaultBase, (*tpl_)[0].Idx);
        ext.Set(0, 0, match_trans_prob * match_emission_prob);
        ext.FinishEditingColumn(0, 0, 1);
    }
}

template <typename Derived>
Recursor<Derived>::Recursor(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr,
                            const double scoreDiff)
    : AbstractRecursor(std::forward<std::unique_ptr<AbstractTemplate>>(tpl), mr, scoreDiff)
    , emissions_{Derived::EncodeRead(read_)}
{
}

template <typename Derived>
size_t Recursor<Derived>::FillAlphaBeta(M& a, M& b) const throw(AlphaBetaMismatch)
{
    if (tpl_->Length() == 0) throw std::runtime_error("template length is 0, invalid state!");

    FillAlpha(M::Null(), a);
    FillBeta(a, b);

    size_t I = read_.Length();
    size_t J = tpl_->Length();
    int flipflops = 0;
    int maxSize = std::max(100, static_cast<int>(0.5 + REBANDING_THRESHOLD * (I + 1) * (J + 1)));

    // if we use too much space, do at least one more round
    // to take advantage of rebanding
    if (a.UsedEntries() >= maxSize || b.UsedEntries() >= maxSize) {
        FillAlpha(b, a);
        FillBeta(a, b);
        FillAlpha(b, a);
        flipflops += 3;
    }

    const double unweight = UndoCounterWeights(read_.Length());
    double alphaV, betaV;
    while (flipflops <= MAX_FLIP_FLOPS) {
        alphaV = std::log(a(I, J)) + a.GetLogProdScales() + unweight;
        betaV = std::log(b(0, 0)) + b.GetLogProdScales() + unweight;

        if (std::abs(1.0 - alphaV / betaV) <= ALPHA_BETA_MISMATCH_TOLERANCE) break;

        if (flipflops % 2 == 0)
            FillAlpha(b, a);
        else
            FillBeta(a, b);

        ++flipflops;
    }

    if (std::abs(1.0 - alphaV / betaV) > ALPHA_BETA_MISMATCH_TOLERANCE || !std::isfinite(betaV))
        throw AlphaBetaMismatch();

    return flipflops;
}

template <typename Derived>
inline Interval Recursor<Derived>::RowRange(size_t j, const M& matrix) const
{
    int beginRow, endRow;
    std::tie(beginRow, endRow) = matrix.UsedRowRange(j);
    int maxRow = beginRow;
    double maxScore = matrix(maxRow, j);
    int i;

    for (i = beginRow + 1; i < endRow; i++) {
        double score = matrix(i, j);

        if (score > maxScore) {
            maxRow = i;
            maxScore = score;
        }
    }

    double thresholdScore = maxScore / scoreDiff_;

    for (i = beginRow; i < maxRow && matrix(i, j) < thresholdScore; i++)
        ;
    beginRow = i;

    for (i = endRow - 1; i >= maxRow && matrix(i, j) < thresholdScore; i--)
        ;
    endRow = i + 1;

    return Interval(beginRow, endRow);
}

template <typename Derived>
inline bool Recursor<Derived>::RangeGuide(size_t j, const M& guide, const M& matrix,
                                          size_t* beginRow, size_t* endRow) const
{
    bool useGuide = !(guide.IsNull() || guide.IsColumnEmpty(j));
    bool useMatrix = !(matrix.IsNull() || matrix.IsColumnEmpty(j));

    if (!useGuide && !useMatrix) {
        return false;
    }

    Interval interval(*beginRow, *endRow);
    if (useGuide) {
        interval = RangeUnion(RowRange(j, guide), interval);
    }

    if (useMatrix) {
        interval = RangeUnion(RowRange(j, matrix), interval);
    }

    std::tie(*beginRow, *endRow) = interval;

    return true;
}

}  // namespace Consensus
}  // namespace PacBio
