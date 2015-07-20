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

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <ConsensusCore/Coverage.hpp>
#include <ConsensusCore/Interval.hpp>

using namespace ConsensusCore;  // NOLINT
using ::testing::ElementsAre;
using ::testing::ElementsAreArray;

#define t(a, b) (Interval((a), (b)))

TEST(CoverageTests, CoverageInWindowTest)
{
    int coverage[10];
    int tStart[] = { 1, 2, 3,  8, 10, 15 };
    int tEnd[]   = { 3, 4, 5, 10, 10, 200};

    CoverageInWindow(6, tStart, 6, tEnd, 0, 10, coverage);
    int expectedCoverage1[] = { 0, 1, 2, 2, 1, 0, 0, 0, 1, 1 };
    ASSERT_THAT(coverage, ElementsAreArray(expectedCoverage1, 10));

    CoverageInWindow(6, tStart, 6, tEnd, 10, 10, coverage);
    int expectedCoverage2[] = { 0, 0, 0, 0, 0, 1, 1, 1, 1, 1 };
    ASSERT_THAT(coverage, ElementsAreArray(expectedCoverage2, 10));
}

TEST(CoverageTests, CoveredIntervalsTest)
{
    // CHECK:
    // CoveredIntervals(0, [1,2,3,8,900,2000], [3,4,5,10,1010,20000], 0, 10000) -> ((0, 10000),)
    // CoveredIntervals(1, [1,2,3,8,900,2000], [3,4,5,10,1010,20000], 0, 10000) -> ((1, 5), (8, 10), (900, 1010), (2000, 10000))  // NOLINT
    // CoveredIntervals(2, [1,2,3,8,900,2000], [3,4,5,10,1010,20000], 0, 10000) -> ((2, 4),)
    // CoveredIntervals(3, [1,2,3,8,900,2000], [3,4,5,10,1010,20000], 0, 10000) -> ()

    int tStart[] = { 1, 2, 3,  8,  900,  2000 };
    int tEnd[]   = { 3, 4, 5, 10, 1010, 20000 };
    ASSERT_THAT(CoveredIntervals(0, 6, tStart, 6, tEnd, 0, 10000), ElementsAre(t(0, 10000)));
    ASSERT_THAT(CoveredIntervals(1, 6, tStart, 6, tEnd, 0, 10000), ElementsAre(t(1, 5), t(8, 10), t(900, 1010), t(2000, 10000)));  // NOLINT
    ASSERT_THAT(CoveredIntervals(2, 6, tStart, 6, tEnd, 0, 10000), ElementsAre(t(2, 4)));
    ASSERT_THAT(CoveredIntervals(3, 6, tStart, 6, tEnd, 0, 10000), ElementsAre());

    ASSERT_THAT(CoveredIntervals(0, 6, tStart, 6, tEnd, 100, 9900), ElementsAre(t(100, 10000)));
    ASSERT_THAT(CoveredIntervals(1, 6, tStart, 6, tEnd, 100, 9900), ElementsAre(t(900, 1010), t(2000, 10000)));  // NOLINT
}


TEST(CoverageTests, CoveredIntervalsTest2)
{
    // Regression test
    int tStart[] = { 48853 };
    int tEnd[]   = { 50687 };
    ASSERT_THAT(CoveredIntervals(1, 1, tStart, 1, tEnd, 50000, 500),
                ElementsAre(t(50000, 50500)));
}
