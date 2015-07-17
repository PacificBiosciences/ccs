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

#include <algorithm>
#include <utility>
#include <string>

#include <ConsensusCore/Types.hpp>
#include <ConsensusCore/Quiver/QuiverConfig.hpp>

namespace ConsensusCore {

    namespace detail {

    /// \brief A base class for recursors, providing some functionality
    ///        based on polymorphic virtual private methods.
    template <typename M, typename E, typename C>
    class RecursorBase
    {
    public:  // Types
        typedef M MatrixType;
        typedef E EvaluatorType;
        typedef C CombinerType;

    public:
        //
        // API methods
        //

        /// \brief Calculate the recursion score by "linking" partial alpha and/or
        ///        beta matrices.
        virtual float LinkAlphaBeta(const E& e,
                                    const M& alpha, int alphaColumn,
                                    const M& beta, int betaColumn,
                                    int absoluteColumn) const = 0;

        /// \brief Fill the alpha and beta matrices.
        /// This routine will fill the alpha and beta matrices, ensuring
        /// that the score computed from the alpha and beta recursions are
        /// identical, refilling back-and-forth if necessary.
        virtual int
        FillAlphaBeta(const E& e, M& alpha, M& beta) const
            throw(AlphaBetaMismatchException);

        /// \brief Reband alpha and beta matrices.
        /// This routine will reband alpha and beta to the convex hull
        /// of the maximum path through each and the inputs for column j.
        virtual bool
        RangeGuide(int j, const M& guide, const M& matrix,
                   int* beginRow, int* endRow) const;

        /// \brief Raw FillAlpha, provided primarily for testing purposes.
        ///        Client code should use FillAlphaBeta.
        virtual void FillAlpha(const E& e, const M& guide, M& alpha) const = 0;

        /// \brief Raw FillBeta, provided primarily for testing purposes.
        ///        Client code should use FillAlphaBeta.
        virtual void FillBeta(const E& e, const M& guide, M& beta) const = 0;

        /// \brief Compute two columns of the alpha matrix starting at columnBegin,
        ///        storing the output in ext.
        virtual void ExtendAlpha(const E& e,
                                 const M& alphaIn, int columnBegin,
                                 M& ext, int numExtColumns = 2) const = 0;


        /// \brief Read out the alignment from the computed alpha matrix.
        const PairwiseAlignment* Alignment(const E& e, const M& alpha) const;


        RecursorBase(int movesAvailable, const BandingOptions& banding);
        virtual ~RecursorBase();

    protected:
        int movesAvailable_;
        BandingOptions bandingOptions_;
    };
}}

#include "RecursorBase-inl.hpp"
