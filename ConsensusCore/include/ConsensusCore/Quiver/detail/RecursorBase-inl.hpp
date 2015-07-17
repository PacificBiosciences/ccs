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

// Author: Lance Hepler

#pragma once

#include <algorithm>
#include <boost/tuple/tuple.hpp>
#include <utility>

#include <ConsensusCore/Interval.hpp>

namespace ConsensusCore {
namespace detail {

    template<typename M>
    inline Interval RowRange(int j, const M& matrix,
                             float scoreDiff)
    {
        int beginRow, endRow;
        boost::tie(beginRow, endRow) = matrix.UsedRowRange(j);
        int maxRow = beginRow;
        float maxScore = matrix(maxRow, j);
        int i;

        for (i = beginRow + 1; i < endRow; i++)
        {
            float score = matrix(i, j);

            if (score > maxScore)
            {
                maxRow = i;
                maxScore = score;
            }
        }

        float thresholdScore = maxScore - scoreDiff;

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

    template<typename M, typename E, typename C>
    inline bool
    RecursorBase<M, E, C>::RangeGuide(int j, const M& guide, const M& matrix,
                                      int* beginRow, int* endRow) const
    {
        bool useGuide = !(guide.IsNull() || guide.IsColumnEmpty(j));
        bool useMatrix = !(matrix.IsNull() || matrix.IsColumnEmpty(j));

        if (!useGuide && !useMatrix)
        {
            return false;
        }

        float scoreDiff = bandingOptions_.ScoreDiff;
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
}}
