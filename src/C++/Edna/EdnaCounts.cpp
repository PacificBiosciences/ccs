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

// Author: Patrick Marks

#include <xmmintrin.h>
#include <pmmintrin.h>
#include <cassert>
#include <cfloat>
#include <iostream>
#include <string>
#include <cmath>
#include <climits>
#include <utility>

#include <ConsensusCore/Edna/EdnaCounts.hpp>
#include <ConsensusCore/Edna/EdnaEvaluator.hpp>
#include <ConsensusCore/Features.hpp>
#include <ConsensusCore/Interval.hpp>
#include <ConsensusCore/Matrix/SparseMatrix.hpp>
#include <ConsensusCore/Quiver/detail/Combiner.hpp>
#include <ConsensusCore/Quiver/detail/RecursorBase.hpp>
#include <ConsensusCore/Quiver/MutationScorer.hpp>
#include <ConsensusCore/Quiver/QuiverConfig.hpp>
#include <ConsensusCore/Quiver/SimpleRecursor.hpp>
#include <ConsensusCore/Types.hpp>
#include <ConsensusCore/Utils.hpp>

using std::min;
using std::max;

#define NEG_INF -FLT_MAX

namespace ConsensusCore
{
    INLINE_CALLEES void
    EdnaCounts::DoCount(Feature<int> channelRead,
                             EdnaEvaluator& eval,
                             MutationScorer<SparseSseEdnaRecursor>& scorer,
                             int j1, int j2, float *results)
    {
        const SparseMatrixF *alpha = scorer.Alpha();
        const SparseMatrixF *beta = scorer.Beta();

        int usedBegin, usedEnd;
        boost::tie(usedBegin, usedEnd) = RangeUnion(alpha->UsedRowRange(j1),
                                                    beta->UsedRowRange(j2));

        for (int k = 0; k < 5; k++)
            results[k] = NEG_INF;

        for (int i = usedBegin; i < usedEnd; i++)
        {
            results[0] = detail::logAdd(results[0],
                                        alpha->Get(i, j1) +
                                        eval.ScoreMove(j1, j2, 0) +
                                        beta->Get(i, j2));
        }

        int nRows = alpha->Rows();
        int usedCap = usedEnd < nRows - 1 ? usedEnd : nRows - 1;

        for (int i = usedBegin; i < usedCap; i++)
        {
            int readBase = channelRead[i];
            results[readBase] = detail::logAdd(results[readBase],
                                               alpha->Get(i, j1) +
                                               eval.ScoreMove(j1, j2, readBase) +
                                               beta->Get(i+1, j2));
        }
    }
}
