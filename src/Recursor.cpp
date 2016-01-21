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

#include <algorithm>
#include <climits>
#include <utility>

#include "Recursor.h"

using std::min;
using std::max;

// TODO(dalexander): put these into a RecursorConfig struct
#define MAX_FLIP_FLOPS 5
#define ALPHA_BETA_MISMATCH_TOLERANCE \
    .001  // TODO: Hmmm... not sure what the heck to do about these...
#define REBANDING_THRESHOLD 0.04

typedef std::pair<size_t, size_t> Interval;

namespace PacBio {
namespace Consensus {
namespace {  // anonymous

inline Interval RangeUnion(const Interval& range1, const Interval& range2)
{
    return Interval(std::min(range1.first, range2.first), std::max(range1.second, range2.second));
}

inline Interval RangeUnion(const Interval& range1, const Interval& range2, const Interval& range3,
                           const Interval& range4)
{
    return RangeUnion(RangeUnion(range1, range2), RangeUnion(range3, range4));
}

double Combine(const double a, const double b) { return a + b; }
}  // namespace anonymous

void Recursor::FillAlpha(const M& guide, M& alpha) const
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
    auto prevTransProbs = TemplatePosition{'-', 0, 0, 0, 0};
    char prevTplBase = 'N';

    for (int j = 1; j < J; ++j)  // Note due to offset with reads and otherwise, this is ugly-ish
    {
        // Load up the transition parameters for this context
        
        auto currTransProbs = (*tpl_)[j - 1];
        auto currTplBase = currTransProbs.Base;
        this->RangeGuide(j, guide, alpha, &hintBeginRow, &hintEndRow);

        size_t requiredEndRow = min(I, hintEndRow);
        size_t i;
        double thresholdScore = 0.0;
        double maxScore = 0.0;
        double score = 0.0;
        alpha.StartEditingColumn(j, hintBeginRow, hintEndRow);

        char nextTplBase = (*tpl_)[j].Base;

        size_t beginRow = hintBeginRow, endRow;
        // Recursively calculate [Probability in last state] * [Probability
        // transition to new state] * [Probability of emission]
        for (i = beginRow; i < I && (score >= thresholdScore || i < requiredEndRow); ++i) {
            const char curReadBase = read_.Seq[i - 1];
            // TODO: Terrible hack right now to emit this guy as teh IQV
            const uint8_t curReadIqv = read_.Seq[i - 1];
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
            double match_prev_and_emission_prob =
                alpha(i - 1, j - 1) *
                tpl_->BaseEmissionPr(MoveType::MATCH, currTplBase, curReadBase);
            if (i == 1 && j == 1) {  // TODO: Remove this branch bottleneck...
                thisMoveScore = match_prev_and_emission_prob;  // Only the emission, since
                                                               // we require a match to
                                                               // start
            } else if (i != 1 && j != 1) {
                thisMoveScore = match_prev_and_emission_prob * prevTransProbs.Match *
                tpl_->CovEmissionPr(MoveType::MATCH, curReadIqv, prevTplBase, currTplBase);
            }
            score = Combine(score, thisMoveScore);

            // Stick or Branch:
            if (i > 1)  // Due to pinning, can't "insert" first or last read base
            {
                const MoveType move =
                    (curReadBase == nextTplBase) ? MoveType::BRANCH : MoveType::STICK;
                double trans_emission_prob =
                    (curReadBase == nextTplBase) ? currTransProbs.Branch : currTransProbs.Stick;
                thisMoveScore = alpha(i - 1, j) * trans_emission_prob *
                                tpl_->BaseEmissionPr(move, nextTplBase, curReadBase) *
                                tpl_->CovEmissionPr(move, curReadIqv, currTplBase, nextTplBase);
                score = Combine(score, thisMoveScore);
            }

            // Deletion:
            if (j > 1)  // Due to pinning, can't "delete" first or last template bp
            {
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
        auto currTplBase = (*tpl_)[J - 1].Base;
        auto prevTplBase = (*tpl_)[J-2].Base;
        // end in the homopolymer state for now.
        auto likelihood = alpha(I - 1, J - 1) *
                          tpl_->BaseEmissionPr(MoveType::MATCH, currTplBase, read_.Seq[I - 1]) *
                          tpl_->CovEmissionPr(MoveType::MATCH, read_.Seq[I - 1], prevTplBase, currTplBase);
        alpha.StartEditingColumn(J, I, I + 1);
        alpha.Set(I, J, likelihood);
        alpha.FinishEditingColumn(J, I, I + 1);
    }
}

void Recursor::FillBeta(const M& guide, M& beta) const
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
        const char nextTplBase = nextTplPos.Base;
        const auto currTransProbs = (*tpl_)[j - 1];

        this->RangeGuide(j, guide, beta, &hintBeginRow, &hintEndRow);

        size_t requiredBeginRow = std::max(hintBeginRow, static_cast<size_t>(0));

        beta.StartEditingColumn(j, hintBeginRow, hintEndRow);

        int i;
        double score = 0.0;
        double thresholdScore = 0.0;
        double maxScore = 0.0;

        int beginRow, endRow = hintEndRow;
        for (i = endRow - 1; i > 0 && (score >= thresholdScore || i >= requiredBeginRow); --i) {
            const char nextReadBase = read_.Seq[i];
            const uint8_t nextReadIqv = read_.Seq[i];
            double thisMoveScore;
            score = 0.0;

            // We are going to pre-catch this because it determine both the
            // match score
            // and whether a stick/branch.
            bool nextBasesMatch = nextReadBase == nextTplBase;

            // Match
            auto match_prev_emission_prob =
                beta(i + 1, j + 1) *
                tpl_->BaseEmissionPr(MoveType::MATCH, nextTplBase, nextReadBase);
            if ((i + 1) < I) {
                thisMoveScore = match_prev_emission_prob * currTransProbs.Match *
                                tpl_->CovEmissionPr(MoveType::MATCH, nextReadIqv, currTransProbs.Base, nextTplBase);
                score = Combine(score, thisMoveScore);
            } else if ((i + 1) == I && ((j + 1) == J)) {
                thisMoveScore = match_prev_emission_prob *
                                tpl_->CovEmissionPr(MoveType::MATCH,
                                                    nextReadIqv, currTransProbs.Base, nextTplBase);  // TODO: Redundant on first pass?
                score = Combine(score, thisMoveScore);
            }

            // Stick or Branch:
            //  can only transition to an insertion for the 2nd to last read
            //  base
            if (0 < i && i < I) {
                const MoveType move = nextBasesMatch ? MoveType::BRANCH : MoveType::STICK;
                auto trans_emission_prob =
                    nextBasesMatch ? currTransProbs.Branch : currTransProbs.Stick;
                thisMoveScore = beta(i + 1, j) * trans_emission_prob *
                                tpl_->BaseEmissionPr(move, nextTplBase, nextReadBase) *
                                tpl_->CovEmissionPr(move, nextReadIqv, currTransProbs.Base, nextTplBase);
                score = Combine(score, thisMoveScore);
            }

            // Deletion:
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
            tpl_->BaseEmissionPr(MoveType::MATCH, (*tpl_)[0].Base, read_.Seq[0]);
        // NO COVARIATE EMISSION FOR FIRST ONE
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
double Recursor::LinkAlphaBeta(const M& alpha, size_t alphaColumn, const M& beta, size_t betaColumn,
                               size_t absoluteColumn) const
{
    const size_t I = read_.Length();

    assert(alphaColumn > 1 && absoluteColumn > 1);
    assert(absoluteColumn <= tpl_->Length());

    int usedBegin, usedEnd;
    std::tie(usedBegin, usedEnd) =
        RangeUnion(alpha.UsedRowRange(alphaColumn - 2), alpha.UsedRowRange(alphaColumn - 1),
                   beta.UsedRowRange(betaColumn), beta.UsedRowRange(betaColumn + 1));

    double v = 0.0, thisMoveScore;

    auto currTplParams = (*tpl_)[absoluteColumn - 1];
    auto currTplBase = currTplParams.Base;

    auto prevTplParams = (*tpl_)[absoluteColumn - 2];

    for (int i = usedBegin; i < usedEnd; i++) {
        if (i < I) {
            const char readBase = read_.Seq[i];
            const uint8_t readIqv = read_.Seq[i];
            double match_prob =
                prevTplParams.Match * tpl_->BaseEmissionPr(MoveType::MATCH, currTplBase, readBase);
            // Incorporate
            thisMoveScore = alpha(i, alphaColumn - 1) * match_prob * beta(i + 1, betaColumn) *
                            tpl_->CovEmissionPr(MoveType::MATCH, readIqv, prevTplParams.Base, currTplParams.Base);
            v = Combine(v, thisMoveScore);
        }

        // Delete:
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
void Recursor::ExtendAlpha(const M& alpha, size_t beginColumn, M& ext, size_t numExtColumns) const
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

    for (size_t extCol = 0; extCol < numExtColumns; extCol++) {
        size_t j = beginColumn + extCol;
        size_t beginRow, endRow;

        //
        // If this extend is contained within the column bounds of
        // the original alpha, we use the row range that was
        // previously determined.  Otherwise start at alpha's last
        // UsedRow beginRow and go to the end.
        //
        // BULLSHIT! If there was a deletion or insertion, the row range for the
        // previous
        // column, not the column of interest will be used.
        // TODO: ERROR! Fix this. Temporary hack is to merge the columns in
        // front and behind.
        // Still totally broken.
        if (j < tpl_->Length()) {
            std::tie(beginRow, endRow) = alpha.UsedRowRange(j);
            size_t pBegin, pEnd, nBegin, nEnd;
            if (j > 0) {
                std::tie(pBegin, pEnd) = alpha.UsedRowRange(j - 1);
                beginRow = std::min(beginRow, pBegin);
                endRow = std::max(endRow, pEnd);
            }
            if ((j + 1) < tpl_->Length()) {
                std::tie(nBegin, nEnd) = alpha.UsedRowRange(j + 1);
                beginRow = std::min(beginRow, nBegin);
                endRow = std::max(endRow, nEnd);
            }
        } else {
            beginRow = alpha.UsedRowRange(alpha.Columns() - 1).first;
            endRow = alpha.Rows();
        }
        ext.StartEditingColumn(extCol, beginRow, endRow);

        int i;
        double score = 0.0;

        // Grab values that will be useful for the whole column
        auto currTplParams = (*tpl_)[j - 1];
        auto currTplBase = currTplParams.Base;
        TemplatePosition prevTplParams{'-', 0, 0, 0, 0};
        if (j > 1) {
            prevTplParams = (*tpl_)[j - 2];
        }
        char nextTplBase;
        if (j != maxLeftMovePossible) {
            nextTplBase = (*tpl_)[j].Base;
        }

        for (i = beginRow; i < endRow; i++) {
            const char currReadBase = read_.Seq[i - 1];
            const uint8_t currReadIqv = read_.Seq[i - 1];
            double thisMoveScore = 0.0;

            // Match:
            if (i > 0 && j > 0) {
                double prev = extCol == 0 ? alpha(i - 1, j - 1) : ext(i - 1, extCol - 1);
                auto emission_prob =
                    tpl_->BaseEmissionPr(MoveType::MATCH, currTplBase, currReadBase);
                
                if (i == 1 && j == 1) {             // TODO: Remove this branch bottleneck...
                    thisMoveScore = emission_prob;  // prev should be 1, so no
                                                    // need for explicit prev +
                                                    // e.Match_Just_Emission
                } else if (i < maxDownMovePossible && j < maxLeftMovePossible) {
                    thisMoveScore = prev * prevTplParams.Match * emission_prob *
                    tpl_->CovEmissionPr(MoveType::MATCH, currReadIqv,
                                        prevTplParams.Base, currTplParams.Base);
                } else if (i == maxDownMovePossible && j == maxLeftMovePossible) {
                    thisMoveScore = prev * emission_prob *
                    tpl_->CovEmissionPr(MoveType::MATCH, currReadIqv,
                                        prevTplParams.Base, currTplParams.Base);;
                }
                score = thisMoveScore;
            }

            // Stick or Branch:
            if (i > 1 && i < maxDownMovePossible && j != maxLeftMovePossible) {
                const MoveType move =
                    (nextTplBase == currReadBase) ? MoveType::BRANCH : MoveType::STICK;
                auto insert_emission_prob =
                    (nextTplBase == currReadBase) ? currTplParams.Branch : currTplParams.Stick;
                thisMoveScore = ext(i - 1, extCol) * insert_emission_prob *
                                tpl_->BaseEmissionPr(move, nextTplBase, currReadBase) *
                                tpl_->CovEmissionPr(move, currReadIqv, currTplBase, nextTplBase);
                score = Combine(score, thisMoveScore);
            }

            // Delete:
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
void Recursor::ExtendBeta(const M& beta, size_t lastColumn, M& ext, int lengthDiff) const
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

    for (int j = lastColumn; j + numExtColumns > lastColumn; j--) {
        /* Convert from old template to new template coordinates.
           lengthDiff will be 0 for substitution, -1 for deletion and +1 for
           insertion
         */
        int jp = j + lengthDiff;
        // What is the current extension column we are adding data into.
        int extCol = lastExtColumn - (lastColumn - j);
        int beginRow, endRow;
        if (j < 0) {
            beginRow = 0;
            endRow = beta.UsedRowRange(0).second;
        } else {
            std::tie(beginRow, endRow) = beta.UsedRowRange(j);
            int pBegin, pEnd, nBegin, nEnd;
            if ((j - 1) >= 0) {
                std::tie(pBegin, pEnd) = beta.UsedRowRange(j - 1);
                beginRow = std::min(beginRow, pBegin);
                endRow = std::max(endRow, pEnd);
            }
            if ((j + 1) < tpl_->Length()) {
                std::tie(nBegin, nEnd) = beta.UsedRowRange(j + 1);
                beginRow = std::min(beginRow, nBegin);
                endRow = std::max(endRow, nEnd);
            }
        }

        ext.StartEditingColumn(extCol, beginRow, endRow);

        // Load up useful values referenced throughout the column.
        auto nextTplParams = (*tpl_)[jp];
        char nextTplBase = nextTplParams.Base;

        TemplatePosition currTplParams{'-', 0, 0, 0, 0};
        if (jp > 0) currTplParams = (*tpl_)[jp - 1];

        for (int i = endRow - 1; i >= beginRow; i--) {
            char nextReadBase = '-';
            unsigned char nextReadIqv = 0;
            if (i < I) {
                nextReadBase = read_.Seq[i];
                nextReadIqv = read_.Seq[i];
            }
            double thisMoveScore = 0.0;
            double score = 0.0;

            bool nextBasesMatch = nextReadBase == nextTplBase;

            // Incorporation:
            // TODO: Remove these checks, we should always be on the left side
            // of the matrix....
            if (i < I && j < J) {
                double next =
                    (extCol == lastExtColumn) ? beta(i + 1, j + 1) : ext(i + 1, extCol + 1);
                double emission_prob =
                    tpl_->BaseEmissionPr(MoveType::MATCH, nextTplBase, nextReadBase);

                // First and last have to start with an emission
                if (((i + 1) == I && (jp + 1) == J) || (i == 0 && j == firstColumn))
                    thisMoveScore = next * emission_prob;
                else if (j > firstColumn && i > 0)
                    thisMoveScore = next * currTplParams.Match * emission_prob *
                    tpl_->CovEmissionPr(MoveType::MATCH, nextReadIqv, currTplParams.Base, nextTplBase);
                score = Combine(score, thisMoveScore);
            }

            // Stick or branch
            if (0 < i && i < I && firstColumn < j) {
                const MoveType move = nextBasesMatch ? MoveType::BRANCH : MoveType::STICK;
                double insert_trans_emission_prob =
                    nextBasesMatch ? currTplParams.Branch : currTplParams.Stick;
                thisMoveScore = ext(i + 1, extCol) * insert_trans_emission_prob *
                                tpl_->BaseEmissionPr(move, nextTplBase, nextReadBase) *
                                tpl_->CovEmissionPr(move, nextReadIqv, currTplParams.Base, nextTplBase);
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
            tpl_->BaseEmissionPr(MoveType::MATCH, (*tpl_)[0].Base, read_.Seq[0]);
        ext.Set(0, 0, match_trans_prob * match_emission_prob);
        ext.FinishEditingColumn(0, 0, 1);
    }
}

Recursor::Recursor(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr,
                   const double scoreDiff)
    : tpl_{std::forward<std::unique_ptr<AbstractTemplate>>(tpl)}
    , read_{mr}
    , scoreDiff_{exp(scoreDiff)}
{
}

size_t Recursor::FillAlphaBeta(M& a, M& b) const throw(AlphaBetaMismatch)
{
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

    double alphaV, betaV;
    while (flipflops <= MAX_FLIP_FLOPS) {
        alphaV = std::log(a(I, J)) + a.GetLogProdScales();
        betaV = std::log(b(0, 0)) + b.GetLogProdScales();

        if (std::abs(1.0 - alphaV / betaV) <= ALPHA_BETA_MISMATCH_TOLERANCE) break;

        if (flipflops % 2 == 0)
            FillAlpha(b, a);
        else
            FillBeta(a, b);

        ++flipflops;
    }

    if (std::abs(1.0 - alphaV / betaV) > ALPHA_BETA_MISMATCH_TOLERANCE) throw AlphaBetaMismatch();

    return flipflops;
}

inline Interval Recursor::RowRange(size_t j, const M& matrix) const
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

inline bool Recursor::RangeGuide(size_t j, const M& guide, const M& matrix, size_t* beginRow,
                                 size_t* endRow) const
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
