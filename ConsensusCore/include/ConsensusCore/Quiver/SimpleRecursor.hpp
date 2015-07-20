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

#include <ConsensusCore/Matrix/DenseMatrix.hpp>
#include <ConsensusCore/Matrix/SparseMatrix.hpp>
#include <ConsensusCore/Quiver/QvEvaluator.hpp>
#include <ConsensusCore/Quiver/detail/Combiner.hpp>
#include <ConsensusCore/Quiver/detail/RecursorBase.hpp>

namespace ConsensusCore {

    template <typename M, typename E, typename C>
    class SimpleRecursor : public detail::RecursorBase<M, E, C>
    {
    public:
        void FillAlpha(const E& e, const M& guide, M& alpha) const;
        void FillBeta(const E& e, const M& guide, M& beta) const;

        float LinkAlphaBeta(const E& e,
                            const M& alpha, int alphaColumn,
                            const M& beta, int betaColumn,
                            int absoluteColumn) const;

        void ExtendAlpha(const E& e,
                         const M& alpha, int beginColumn,
                         M& ext, int numExtColumns = 2) const;

        void ExtendBeta(const E& e,
                        const M& beta, int endColumn,
                        M& ext, int numExtColumns = 2,
                        int lengthDiff = 0) const;


    public:
        //
        // Constructors
        //
        SimpleRecursor(int movesAvailable, const BandingOptions& banding);
    };

    typedef SimpleRecursor<DenseMatrixF,
                           QvEvaluator,
                           detail::ViterbiCombiner> SimpleQvRecursor;

    typedef SimpleRecursor<DenseMatrixF,
                           QvEvaluator,
                           detail::SumProductCombiner> SimpleQvSumProductRecursor;

    typedef SimpleRecursor<SparseMatrixF,
                           QvEvaluator,
                           detail::ViterbiCombiner> SparseSimpleQvRecursor;

    typedef SimpleRecursor<SparseMatrixF,
                           QvEvaluator,
                           detail::SumProductCombiner> SparseSimpleQvSumProductRecursor;
}
