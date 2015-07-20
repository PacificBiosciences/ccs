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

#pragma once

#include <ConsensusCore/Logging.hpp>
#include <ConsensusCore/Read.hpp>
#include <ConsensusCore/Matrix/ScaledMatrix.hpp>
#include <ConsensusCore/Arrow/detail/Combiner.hpp>
#include <ConsensusCore/Arrow/ArrowConfig.hpp>
#include <ConsensusCore/Arrow/TemplateParameterPair.hpp>

namespace ConsensusCore {
namespace Arrow {

    template <typename M, typename C>
    class SimpleRecursor
    {

    public:
        
        // Types
        typedef M MatrixType;
        typedef C CombinerType;
        typedef ModelParams ParamsType;
        
        ArrowRead read_;
        WrappedTemplateParameterPair tpl_;
        ModelParams params_;
        
        
        /// \brief Fill the alpha and beta matrices.
        /// This routine will fill the alpha and beta matrices, ensuring
        /// that the score computed from the alpha and beta recursions are
        /// identical, refilling back-and-forth if necessary.
        int
        FillAlphaBeta( M& alpha, M& beta) const
        throw(AlphaBetaMismatchException);
        
        /**
         Fill in the alpha matrix.  This matrix has the read run along the rows, and the
         template run along the columns.  The first row and column do not correspond to
         a template position.  Therefore the match represented at position (i,j) corresponds
         to a match between template positions (i+1, j+1).
         
         The alpha matrix is the "Forward" matrix used in the forward/backward algorithm.
         The i,j position of the matrix represents the probability of all paths up
         to the point where the ith read position and jth template have been "emitted."
         The matrix is calculated recursively by examining all possible transitions
         into (i,j), and calculating the probability we were in the previous state,
         times the probability of a transition into (i,j) times the probability of
         emitting the observation that corresponds to (i,j). All probabilities are
         calculated and stored as LOG values.
         
         Note that in doing this calculation, in order to work with di-nucleotide contexts, we
         require that the first and last transition be a match.  In other words the start and end of
         the read and template are "pinned" to each other.
         
         //TODO: Verify memory is initialized to 0!
         
         @param guide An object that helps inform how to select the size of "bands" for the
         banded algorithm used.  This is typically the beta matrix if we are "repopulating" the matrix.
         @param alpha The matrix to be filled.
         */

        void FillAlpha(const M& guide, M& alpha) const;
        
        /**
         Fill the Beta matrix, the backwards half of the forward-backward algorithm.
         This represents the probability that starting from the (i,j) state, the combined
         probability of transitioning out and following all paths through to the end.
         That is, we need to calculate transition from state and emit from next state for each
         
         In combination with the Alpha matrix, this allows us to calculate all paths that
         pass through the (i,j) element, as exp(Alpha(i,j) + Beta(i,j))
         
         All probabilities stored in the matrix are stored as NON-LOGGED probabilities.
         
         @param e The evaluator, such as QvEvaluator
         @param M the guide matrix for banding (this needs more documentation)
         @param beta The Beta matrix, stored as either a DenseMatrix or a SparseMatrix.
         */

        void FillBeta(const M& guide, M& beta) const;
        
        /// \brief Calculate the recursion score by "linking" partial alpha and/or
        ///        beta matrices.
        double LinkAlphaBeta(const M& alpha, int alphaColumn,
                            const M& beta, int betaColumn,
                            int absoluteColumn) const;

        void ExtendAlpha(const M& alpha, int beginColumn,
                         M& ext, int numExtColumns = 2) const;

        void ExtendBeta(const M& beta, int endColumn,
                        M& ext, int lengthDiff = 0) const;


        //
        // Constructors
        //
        SimpleRecursor(const ModelParams& params,
                       const ArrowRead& read,
                       const WrappedTemplateParameterPair& wtpp,
                       const BandingOptions& banding);
        
#pragma mark Row guide functions for banding optimizations.
    private:
        Interval RowRange(int j, const M& matrix, double scoreDiff) const;
        
        /// \brief Reband alpha and beta matrices.
        /// This routine will reband alpha and beta to the convex hull
        /// of the maximum path through each and the inputs for column j.
        bool RangeGuide(int j, const M& guide, const M& matrix, int* beginRow, int* endRow) const;
        BandingOptions bandingOptions_;
    };

    typedef SimpleRecursor<ConsensusCore::ScaledSparseMatrixD,
                           detail::SumProductCombiner> ArrowRecursor;
}
}
