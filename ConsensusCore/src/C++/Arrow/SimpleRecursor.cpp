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
#include <boost/tuple/tuple.hpp>
#include <climits>
#include <utility>

#include <ConsensusCore/Matrix/ScaledMatrix.hpp>
#include <ConsensusCore/Arrow/detail/Combiner.hpp>
#include <ConsensusCore/Arrow/SimpleRecursor.hpp>
#include <ConsensusCore/Interval.hpp>
#include <ConsensusCore/Utils.hpp>

using std::min;
using std::max;

// TODO(dalexander): put these into a RecursorConfig struct
#define MAX_FLIP_FLOPS                  5
#define ALPHA_BETA_MISMATCH_TOLERANCE   .001 // TODO: Hmmm... not sure what the heck to do about these...
#define REBANDING_THRESHOLD             0.04


namespace ConsensusCore {
namespace Arrow {

    template<typename M, typename C>
    void
    SimpleRecursor<M, C>::FillAlpha(const M& guide, M& alpha) const
    {
        // We are pinning, so should never go all the way to the end of the read/template
        // But our matrix indexing is one off the model/outcome indexing
        // so the match in (1,1) corresponds to a pairing between Model[0]/Outcome[0]
        int I = read_.Length();
        int J = tpl_.Length();

        assert(alpha.Rows() == I + 1 && alpha.Columns() == J + 1);
        assert(guide.IsNull() ||
               (guide.Rows() == alpha.Rows() && guide.Columns() == alpha.Columns()));


        // Initial condition, we always start with a match
        alpha.StartEditingColumn(0, 0, 1);
        alpha.Set(0, 0, 1.0);
        alpha.FinishEditingColumn(0, 0, 1);
        // End initial conditions


        // NOPE:
        // int hintBeginRow = 1, hintEndRow = I - 1;
        int hintBeginRow = 1, hintEndRow = 1;
        auto prevTransProbs = TransitionParameters();
        const ModelParams& params = params_;
        auto score_diff_natural_scale = std::exp(this->bandingOptions_.ScoreDiff);


        for (int j = 1; j < J; ++j) // Note due to offset with reads and otherwise, this is ugly-ish
        {
            // Load up the transition parameters for this context
            auto currentTempPos =  tpl_.GetTemplatePosition(j - 1);
            auto curTransProbs = currentTempPos.second;
            auto curTempBase = currentTempPos.first;

            this->RangeGuide(j, guide, alpha, &hintBeginRow, &hintEndRow);

            int requiredEndRow = min(I , hintEndRow);
            int i;
            double thresholdScore = 0.0;
            double maxScore = 0.0;
            double score = 0.0;
            alpha.StartEditingColumn(j, hintBeginRow, hintEndRow);

            char nextTplBase = tpl_.GetTemplatePosition(j).first;

            int beginRow = hintBeginRow, endRow;
            // Recursively calculate [Probability in last state] * [Probability transition to new state] * [Probability of emission]
            for (i = beginRow;
                 i < I && (score >= thresholdScore || i < requiredEndRow);
                 ++i)
            {
                auto curReadBase = read_.Features[i - 1]; //
                unsigned char curReadIqv = read_.Features.InsQv[i - 1];
                double thisMoveScore = 0.0;
                score = 0.0;
                // Match:
                    /* Important!  Note that because we require the initial state to be a match,
                       when i = 1 and j = 1 the match transition probability must be 1, since no other options
                       are allowed.  Similarly, the probability for the match probability to the end base should be 1.

                       Note that for the first "match" between a read and template, we have no choice but to
                       hard code it to 1, as there is no defined transition probability for a dinucleotide context.

                      ***********  EDGE_CONDITION ************
                     */
                    double match_prev_and_emission_prob = alpha(i - 1, j - 1) * (curReadBase == curTempBase ? params.PrNotMiscall : params.PrThirdOfMiscall);
                    if (i == 1 && j == 1) { //TODO: Remove this branch bottleneck...
                        thisMoveScore =  match_prev_and_emission_prob; // Only the emission, since we require a match to start
                    }
                    else if (i != 1 && j != 1) {
                        thisMoveScore = match_prev_and_emission_prob * prevTransProbs.Match;
                    }
                    score = C::Combine(score, thisMoveScore * params.MatchIqvPmf[curReadIqv]);

                // Stick or Branch:
                    if (i > 1) // Due to pinning, can't "insert" first or last read base
                    {
                        auto trans_emission_prob = curReadBase == nextTplBase ? curTransProbs.Branch : (curTransProbs.Stick / 3.0);
                        thisMoveScore = alpha(i - 1, j) * trans_emission_prob * params.InsertIqvPmf[curReadIqv];
                        score = C::Combine(score, thisMoveScore);
                    }

                // Deletion:
                    if (j > 1) // Due to pinning, can't "delete" first or last template bp
                    {
                        thisMoveScore = alpha(i, j - 1) * prevTransProbs.Deletion;
                        score = C::Combine(score, thisMoveScore);
                    }

                //  Save score
                    alpha.Set(i, j, score);
                    if (score > maxScore)
                    {
                        maxScore = score;
                        thresholdScore = maxScore / score_diff_natural_scale;
                    }
            }
            endRow = i;
            alpha.FinishEditingColumn(j, beginRow, endRow);
            prevTransProbs = curTransProbs;
            // Now, revise the hints to tell the caller where the mass of the
            // distribution really lived in this column.
            hintEndRow = endRow;
            for (i = beginRow; i < endRow && alpha(i, j) < thresholdScore; ++i);
            hintBeginRow = i;
        }
        auto currentTempPos =  tpl_.GetTemplatePosition(J-1);
        auto curTempBase = currentTempPos.first;

        /* Now fill out the probability in the last pinned position.
         * We require that we end in a match.
         * search for the term EDGE_CONDITION to find a comment with more information */
        auto match_emission =  (read_.Features[I-1] == curTempBase ? params.PrNotMiscall : params.PrThirdOfMiscall);
        auto likelihood = alpha(I - 1, J - 1) * match_emission * params.MatchIqvPmf[read_.Features.InsQv[I-1]];
        alpha.StartEditingColumn(J, I, I + 1);
        alpha.Set(I, J, likelihood);
        alpha.FinishEditingColumn(J, I, I + 1);
        //DumpMatrix(alpha);
    }

    template<typename M, typename C>
    void
    SimpleRecursor<M, C>::FillBeta(const M& guide, M& beta) const
    {
        int I = read_.Length();
        int J = tpl_.Length();

        assert(beta.Rows() == I + 1 && beta.Columns() == J + 1);
        assert(guide.IsNull() || (guide.Rows() == beta.Rows() && guide.Columns() == beta.Columns()));

        //Setup initial condition, at the end we are one
        beta.StartEditingColumn(J, I, I + 1);
        beta.Set(I, J, 1.0);
        beta.FinishEditingColumn(J, I, I + 1);

        const ModelParams& params = params_;
        auto score_diff_natural_scale = std::exp(this->bandingOptions_.ScoreDiff);

        // Totally arbitray decision here...
        int hintBeginRow = I, hintEndRow = I;
        // NOPE:
        // int hintBeginRow = std::min(0, I - 45), hintEndRow = I;

        // Recursively calculate [Probability transition to next state] * [Probability of emission at that state] * [Probability from that state]
        for (int j = (J-1); j > 0; --j)
        {
            auto nextTempPos = tpl_.GetTemplatePosition(j);
            auto nextTempBase = nextTempPos.first;

            auto curTempPos = tpl_.GetTemplatePosition(j - 1);
            auto curTransProbs = curTempPos.second;

            this->RangeGuide(j, guide, beta, &hintBeginRow, &hintEndRow);

            int requiredBeginRow = max(0, hintBeginRow);

            beta.StartEditingColumn(j, hintBeginRow, hintEndRow);

            int i;
            double score = 0.0;
            double thresholdScore = 0.0;
            double maxScore = 0.0;

            int beginRow, endRow = hintEndRow;
            for (i = endRow - 1;
                 i > 0 && (score >= thresholdScore || i >= requiredBeginRow);
                 --i)
            {
                auto nextReadBase = read_.Features[i];
                unsigned char nextReadIqv = read_.Features.InsQv[i];
                double thisMoveScore;
                score = 0.0;

                // We are going to pre-catch this because it determine both the match score
                // and whether a stick/branch.
                bool nextBasesMatch = nextReadBase == nextTempBase;


                // Match
                    auto match_prev_emission_prob = beta(i+1, j+1) * (nextBasesMatch ? params.PrNotMiscall : params.PrThirdOfMiscall);
                    if (i < (I-1)) {
                         score = C::Combine(score, match_prev_emission_prob * curTransProbs.Match * params.MatchIqvPmf[nextReadIqv]);
                    }
                    else if (i == (I-1) && (j == (J-1))) {

                        score = C::Combine(score, match_prev_emission_prob * params.MatchIqvPmf[nextReadIqv]); // TODO: Redundant on first pass?
                    }

                // Stick or Branch:
                    if (i < (I - 1) && i > 0) // Can only transition to an insertion for the 2nd to last read base
                    {
                        auto trans_emission_prob = nextBasesMatch ? curTransProbs.Branch : (curTransProbs.Stick / 3.0);
                        thisMoveScore = beta(i + 1, j) * trans_emission_prob * params.InsertIqvPmf[nextReadIqv];
                        score = C::Combine(score, thisMoveScore);
                    }

                // Deletion:
                    if (j < (J - 1) && j > 0)
                    {
                        thisMoveScore = beta(i, j + 1) * curTransProbs.Deletion;
                        score = C::Combine(score, thisMoveScore);
                    }

                // Save score
                beta.Set(i, j, score);

                if (score > maxScore)
                {
                    maxScore = score;
                    thresholdScore = maxScore / score_diff_natural_scale;
                }
            }

            beginRow = i + 1;
            beta.FinishEditingColumn(j, beginRow, endRow);
            //DumpBetaMatrix(beta);
            // Now, revise the hints to tell the caller where the mass of the
            // distribution really lived in this column.
            hintBeginRow = beginRow;
            for (i = endRow;
                 i > beginRow && beta(i - 1, j) < thresholdScore;
                 --i);
            hintEndRow = i;
        }

        beta.StartEditingColumn(0, 0, 1);
        /* Now to fill the top row which must be a match
         * search for the term EDGE_CONDITION to find a comment with more information */

        auto match_emission_prob = (tpl_.GetTemplatePosition(0).first == read_.Features[0]) ? params.PrNotMiscall : params.PrThirdOfMiscall;
        auto match_iqv_emisson_prob = params.MatchIqvPmf[read_.Features.InsQv[0]];
        beta.Set(0, 0, match_emission_prob * beta(1, 1) * match_iqv_emisson_prob);
        beta.FinishEditingColumn(0, 0, 1);
    }

    /// Calculate the recursion score by "stitching" together partial
    /// alpha and beta matrices.  alphaColumn, betaColumn, and
    /// absoluteColumn all refer to the same logical position in the
    /// template, but may have different values if, for instance,
    /// alpha here is a sub-range of the columns of the full alpha
    /// matrix.  Columns betaColumn and betaColumn + 1 of beta will be
    /// read; columns alphaColumn - 1 and alphaColumn - 2 of alpha
    /// will be read.
    template<typename M, typename C>
    double
    SimpleRecursor<M, C>::LinkAlphaBeta(const M& alpha, int alphaColumn,
                                        const M& beta, int betaColumn,
                                        int absoluteColumn) const
    {
        const int I = read_.Length();

        assert(alphaColumn > 1 && absoluteColumn > 1);
        assert(absoluteColumn <= tpl_.Length());

        int usedBegin, usedEnd;
        boost::tie(usedBegin, usedEnd) = \
            RangeUnion(alpha.UsedRowRange(alphaColumn - 2),
                       alpha.UsedRowRange(alphaColumn - 1),
                       beta.UsedRowRange(betaColumn),
                       beta.UsedRowRange(betaColumn + 1));

        double v = 0.0, thisMoveScore;

        auto currentTplPos = tpl_.GetTemplatePosition(absoluteColumn - 1);
        auto currentTplBase = currentTplPos.first;

        auto prevTplPos = tpl_.GetTemplatePosition(absoluteColumn - 2);
        TransitionParameters prevTplParams = prevTplPos.second;

        for (int i = usedBegin; i < usedEnd; i++)
        {
            if (i < I)
            {
                char readBase = read_.Features[i];
                unsigned char readIqv = read_.Features.InsQv[i];
                double match_prob = prevTplParams.Match * (readBase == currentTplBase ? params_.PrNotMiscall : params_.PrThirdOfMiscall);
                // Incorporate
                thisMoveScore = alpha(i, alphaColumn - 1) *
                                match_prob *
                                beta(i + 1, betaColumn) *
                                params_.MatchIqvPmf[readIqv];
                v = C::Combine(v, thisMoveScore);
            }

            // Delete:
            thisMoveScore = alpha(i, alphaColumn - 1) *
                            prevTplParams.Deletion *
                            beta(i, betaColumn);
            v = C::Combine(v, thisMoveScore);
        }

        return ( std::log(v)
               + alpha.GetLogProdScales(0, alphaColumn)
               + beta.GetLogProdScales(betaColumn, beta.Columns()) );
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
    template<typename M, typename C>
    void
    SimpleRecursor<M, C>::ExtendAlpha(const M& alpha, int beginColumn,
                                         M& ext, int numExtColumns) const
    {
        assert(numExtColumns >= 2); // We have to fill at least one
        assert(alpha.Rows() == read_.Length() + 1 &&
               ext.Rows() == read_.Length() + 1); // The read never mutates

        // The new template may not be the same length as the old template.
        // Just make sure that we have anough room to fill out the extend buffer
        assert(beginColumn + 1 < tpl_.Length() + 1);
        assert(ext.Columns() >= numExtColumns);
        assert(beginColumn >= 2);
        // Due to pinning at the end, moves are only possible if less than these positions.
        int maxLeftMovePossible = tpl_.Length();
        int maxDownMovePossible = read_.Length();

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
            // BULLSHIT! If there was a deletion or insertion, the row range for the previous
            // column, not the column of interest will be used.
            // TODO: ERROR! Fix this. Temporary hack is to merge the columns in front and behind.
            // Still totally broken.
            if (j < tpl_.Length())
            {
                boost::tie(beginRow, endRow) = alpha.UsedRowRange(j);
                int pBegin, pEnd, nBegin, nEnd;
                if ( (j-1) >= 0) {
                    boost::tie(pBegin, pEnd) = alpha.UsedRowRange(j-1);
                    beginRow = std::min(beginRow, pBegin);
                    endRow = std::max(endRow, pEnd);
                }
                if ( (j+1) < tpl_.Length() )  {
                    boost::tie(nBegin, nEnd) = alpha.UsedRowRange(j+1);
                    beginRow = std::min(beginRow, nBegin);
                    endRow = std::max(endRow, nEnd);
                }
            }
            else
            {
                beginRow = alpha.UsedRowRange(alpha.Columns() - 1).Begin;
                endRow = alpha.Rows();
            }
            ext.StartEditingColumn(extCol, beginRow, endRow);

            int i;
            double score;

            // Grab values that will be useful for the whole column
            auto currentTplPos = tpl_.GetTemplatePosition(j - 1);
            auto curTplBase = currentTplPos.first;
            auto curTplParams = currentTplPos.second;
            TransitionParameters prevTplParams;
            if (j > 1) {
                prevTplParams = tpl_.GetTemplatePosition(j-2).second;
            }
            char nextTplBase;
            if ( j != maxLeftMovePossible ) {
                nextTplBase = tpl_.GetTemplatePosition(j).first;
            }


            for (i = beginRow; i < endRow; i++)
            {
                auto curReadBase = read_.Features[i - 1];
                unsigned char curReadIqv = read_.Features.InsQv[i - 1];
                double thisMoveScore = 0.0;

                // Match:
                if (i > 0 && j > 0) {
                    double prev = extCol == 0 ? alpha(i - 1, j - 1) : ext(i - 1, extCol - 1);
                    auto emission_prob = curReadBase == curTplBase ? params_.PrNotMiscall : params_.PrThirdOfMiscall;
                    if (i == 1 && j == 1) { //TODO: Remove this branch bottleneck...
                        thisMoveScore = emission_prob; // prev should be 1, so no need for explicit prev + e.Match_Just_Emission
                    }
                    else if (i < maxDownMovePossible && j < maxLeftMovePossible) {
                        thisMoveScore = prev * prevTplParams.Match * emission_prob;
                    }
                    else if (i == maxDownMovePossible && j == maxLeftMovePossible) {
                        thisMoveScore = prev * emission_prob;
                    }
                    score = thisMoveScore * params_.MatchIqvPmf[curReadIqv];
                }

                // Stick or Branch:
                if (i > 1 && i < maxDownMovePossible && j != maxLeftMovePossible)
                {
                    auto insert_emission_prob = (nextTplBase == curReadBase) ? curTplParams.Branch : (curTplParams.Stick / 3.0);
                    thisMoveScore = ext(i - 1, extCol) * insert_emission_prob * params_.InsertIqvPmf[curReadIqv];
                    score = C::Combine(score, thisMoveScore);
                }

                // Delete:
                if (j > 1 && j < maxLeftMovePossible && i != maxDownMovePossible)
                {
                    double prev = extCol == 0 ? alpha(i, j - 1) : ext(i, extCol - 1);
                    thisMoveScore = prev * prevTplParams.Deletion;
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

    /* Note this is a very confusing routine in order to avoid recomputing and
       additional memory allocations.  This routine tries to stick on a beta matrix
       to the original and back trace to the 0,0 position of this extension matrix.
       matrix from the original.  Note that the original beta
       matrix is indexed by the original template positions, while the template
       bases and parameters are now indexed according the the "virtual" template
       to which mutations have been applied.
     */
    // @param lastColumn - Where we
    template<typename M, typename C>
    void
    SimpleRecursor<M, C>::ExtendBeta(const M& beta, int lastColumn,
                                        M& ext, int lengthDiff) const
    {

        int I = read_.Length();
        int J = tpl_.Length();
        // How far back do we have to go until we are at the zero (first) column?
        // we always go all the way back.
        int numExtColumns = lengthDiff + lastColumn + 1;
        int firstColumn = 0 - lengthDiff;
        int lastExtColumn = numExtColumns - 1;

        // The new template may not be the same length as the old template.
        // Just make sure that we have enough room to fill out the extend buffer
        assert(lastColumn + 2 <= J);
        assert(lastColumn < 4); // Since we are only testing mutations of size 1, and the check prior for a beginning mutation is < 3, max = 2 + 1 = 3
        assert(lastColumn >= 0);
        assert(ext.Columns() >= numExtColumns);
        assert(beta.Rows() == I + 1 && ext.Rows() == I + 1);
        assert(std::abs(lengthDiff) < 2);

        for (int j = lastColumn; j > lastColumn - numExtColumns; j--)
        {
            /* Convert from old template to new template coordinates.
               lengthDiff will be 0 for substitution, -1 for deletion and +1 for insertion
             */
            int jp = j + lengthDiff;
            // What is the current extension column we are adding data into.
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
                int pBegin, pEnd, nBegin, nEnd;
                if ( (j-1) >= 0) {
                    boost::tie(pBegin, pEnd) = beta.UsedRowRange(j-1);
                    beginRow = std::min(beginRow, pBegin);
                    endRow = std::max(endRow, pEnd);
                }
                if ( (j+1) < tpl_.Length() ) {
                    boost::tie(nBegin, nEnd) = beta.UsedRowRange(j+1);
                    beginRow = std::min(beginRow, nBegin);
                    endRow = std::max(endRow, nEnd);
                }
            }

            ext.StartEditingColumn(extCol, beginRow, endRow);

            int i;
            double score;


            // Load up useful values referenced throughout the column.
            char nextTplBase;
            auto nextTplPos = tpl_.GetTemplatePosition(jp);
            nextTplBase = nextTplPos.first;

            TransitionParameters curTransParams;
            if (jp > 0) {
                auto curTransPos = tpl_.GetTemplatePosition(jp - 1);
                curTransParams = curTransPos.second;
            }

            for (i = endRow - 1; i >= beginRow; i--)
            {
                char nextReadBase = 'N';
                unsigned char nextReadIqv = 0;
                if (i < I) {
                    nextReadBase = read_.Features[i];
                    nextReadIqv = read_.Features.InsQv[i];
                }
                double thisMoveScore = 0.0;
                score = 0.0;

                bool nextBasesMatch = nextReadBase == nextTplBase;

                // Incorporation:
                // TODO: Remove these checks, we should always be on the left side of the matrix....
                if (i < I && j < J) // Should be I-1, J-1 or removed entirely....
                {
                    double next = (extCol == lastExtColumn) ? beta(i + 1, j + 1) : ext(i + 1, extCol + 1);
                    double emission_prob = nextBasesMatch ? params_.PrNotMiscall : params_.PrThirdOfMiscall;
                    // First and last have to start with an emission
                    if ((i == (I-1) && jp == (J-1)) || (i==0 && j == firstColumn)) {
                        thisMoveScore = next * emission_prob;
                    }
                    else if ( j > firstColumn && i > 0) {
                        thisMoveScore = next * curTransParams.Match * emission_prob;
                    }
                    score = C::Combine(score, thisMoveScore * params_.MatchIqvPmf[nextReadIqv]);
                }

                // Stick or branch
                if (i < (I-1) && i > 0 && j > firstColumn)
                {
                    double insert_trans_emission_prob = nextBasesMatch ? curTransParams.Branch : (curTransParams.Stick / 3.0);
                    thisMoveScore = ext(i + 1, extCol) * insert_trans_emission_prob * params_.InsertIqvPmf[nextReadIqv];
                    score = C::Combine(score, thisMoveScore);
                }

                // Deletion
                if (j < (J-1) && j > firstColumn && i > 0)
                {
                    double next = (extCol == lastExtColumn) ? beta(i, j + 1) : ext(i, extCol + 1);
                    thisMoveScore = next * curTransParams.Deletion;
                    score = C::Combine(score, thisMoveScore);
                }

                ext.Set(i, extCol, score);
            }
            ext.FinishEditingColumn(extCol, beginRow, endRow);
        }
    }


    template<typename M, typename C>
    SimpleRecursor<M, C>::SimpleRecursor(const ModelParams& params,
                                         const ArrowRead& read,
                                         const WrappedTemplateParameterPair& wtpp,
                                         const BandingOptions& banding)
        : read_(read)
        , tpl_(wtpp)
        , params_(params)
        , bandingOptions_(banding)
    {}

    template<typename M, typename C>
    int
    SimpleRecursor<M, C>::FillAlphaBeta(M& a, M& b) const
    throw(AlphaBetaMismatchException)
    {
        FillAlpha(M::Null(), a);
        FillBeta( a, b);

        int I = read_.Length();
        int J = tpl_.Length();
        int flipflops = 0;
        int maxSize = static_cast<int>(0.5 + REBANDING_THRESHOLD * (I + 1) * (J + 1));

        // if we use too much space, do at least one more round
        // to take advantage of rebanding
        if (a.UsedEntries() >= maxSize ||
            b.UsedEntries() >= maxSize)
        {
            FillAlpha(b, a);
            FillBeta(a, b);
            FillAlpha(b, a);
            flipflops += 3;
        }
        double alphaV = std::log(a(I, J)) + a.GetLogProdScales();
        double betaV  = std::log(b(0, 0)) + b.GetLogProdScales();
        while (fabs(alphaV - betaV) > ALPHA_BETA_MISMATCH_TOLERANCE
               && flipflops <= MAX_FLIP_FLOPS)
        {
            if (flipflops % 2 == 0)
            {
                FillAlpha(b, a);
            }
            else
            {
                FillBeta(a, b);
            }
            flipflops++;
        }
        alphaV = std::log(a(I, J)) + a.GetLogProdScales();
        betaV  = std::log(b(0, 0)) + b.GetLogProdScales();
        auto mismatch_percentage = fabs(1.0 - alphaV/betaV);
        if (mismatch_percentage > ALPHA_BETA_MISMATCH_TOLERANCE)
        {
            LDEBUG << "Could not mate alpha, beta.  Read: " << read_.Name << " Tpl: Was wrapped, improve debugging to pring";
            // std::cout << "Mismatch : " << mismatch_percentage << "% Alpha " << alphaV << " Beta " << betaV << std::endl;
            throw AlphaBetaMismatchException();
        }

        return flipflops;
    }

#pragma mark Row guide functions for banding optimizations.
    template<typename M, typename C>
    inline Interval SimpleRecursor<M, C>::RowRange(int j, const M& matrix,  double scoreDiff) const
    {
        int beginRow, endRow;
        boost::tie(beginRow, endRow) = matrix.UsedRowRange(j);
        int maxRow = beginRow;
        double maxScore = matrix(maxRow, j);
        int i;

        for (i = beginRow + 1; i < endRow; i++)
        {
            double score = matrix(i, j);

            if (score > maxScore)
            {
                maxRow = i;
                maxScore = score;
            }
        }

        double thresholdScore = maxScore - scoreDiff;

        for (i = beginRow;
             i < maxRow && matrix(i, j) < thresholdScore;
             i++);
        beginRow = i;

        for (i = endRow - 1;
             i >= maxRow && matrix(i, j) < thresholdScore;
             i--);
        endRow = i + 1;

        return Interval(beginRow, endRow);
    }

    template<typename M, typename C>
    inline bool
    SimpleRecursor<M, C>::RangeGuide(int j, const M& guide, const M& matrix,
                                       int* beginRow, int* endRow) const
    {
        bool useGuide = !(guide.IsNull() || guide.IsColumnEmpty(j));
        bool useMatrix = !(matrix.IsNull() || matrix.IsColumnEmpty(j));

        if (!useGuide && !useMatrix)
        {
            return false;
        }

        double scoreDiff = bandingOptions_.ScoreDiff;
        Interval interval(*beginRow, *endRow);

        if (useGuide)
        {
            interval = RangeUnion(RowRange(j, guide, scoreDiff), interval);
        }

        if (useMatrix)
        {
            interval = RangeUnion(RowRange(j, matrix, scoreDiff), interval);
        }

        boost::tie(*beginRow, *endRow) = interval;

        return true;
    }

    template class SimpleRecursor<ScaledSparseMatrixD, detail::SumProductCombiner>;
}
}
