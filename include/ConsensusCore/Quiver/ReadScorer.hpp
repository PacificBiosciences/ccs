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

#include <ConsensusCore/Read.hpp>
#include <ConsensusCore/Types.hpp>
#include <ConsensusCore/Matrix/SparseMatrix.hpp>
#include <ConsensusCore/Quiver/QuiverConfig.hpp>

#include <string>

namespace ConsensusCore
{
    // A ReadScorer provides a simple but inefficient means of scoring
    // a single read against a template.  This is only to be used for
    // testing, not in production code.

    class ReadScorer
    {
    public:
        explicit ReadScorer(QuiverConfig& config);

        float Score
        (const std::string& tpl, const QvRead& read) const
            throw(AlphaBetaMismatchException);

        const PairwiseAlignment* Align
        (const std::string& tpl, const QvRead& read) const
            throw(AlphaBetaMismatchException);

        const SparseMatrixF* Alpha
        (const std::string& tpl, const QvRead& read) const
            throw(AlphaBetaMismatchException);

        const SparseMatrixF* Beta
        (const std::string& tpl, const QvRead& read) const
            throw(AlphaBetaMismatchException);


    private:
        QuiverConfig _quiverConfig;
    };
}
