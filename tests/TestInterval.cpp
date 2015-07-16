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

#include <iostream>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <pacbio/ccs/Interval.h>
#include <pacbio/ccs/IntervalTree.h>

using namespace PacBio::CCS;

TEST(IntervalTest, Merging)
{
    IntervalTree tree;

    tree.Insert(Interval(1, 3));
    tree.Insert(Interval(3, 5));

    EXPECT_EQ(tree.size(), 1);

    for (const auto& i : tree)
    {
        EXPECT_EQ(i.Right(), 5);
        EXPECT_EQ(i.Left(),  1);
    }
}

TEST(IntervalTest, Merging2)
{
    IntervalTree tree;

    tree.Insert(Interval(1, 3));
    tree.Insert(Interval(5, 7));
    tree.Insert(Interval(9, 11));

    EXPECT_EQ(tree.size(), 3);

    tree.Insert(Interval(3, 9));

    EXPECT_EQ(tree.size(), 1);

    for (const auto& i : tree)
    {
        EXPECT_EQ(i.Right(), 11);
        EXPECT_EQ(i.Left(),  1);
    }
}

TEST(IntervalTest, Merging3)
{
    IntervalTree tree;

    tree.Insert(Interval(1, 3));
    tree.Insert(Interval(5, 6));
    tree.Insert(Interval(4, 6));

    EXPECT_EQ(tree.size(), 2);
}

TEST(IntervalTest, Iteration)
{
    Interval interval(0, 11);
    size_t i = interval.Left();

    for (const auto j : interval)
    {
        EXPECT_EQ(j, i++);
    }
}

TEST(IntervalTest, Gaps)
{
    IntervalTree tree;

    tree.Insert(Interval(1, 3));
    tree.Insert(Interval(5, 7));
    tree.Insert(Interval(9, 11));

    IntervalTree gaps = tree.Gaps();

    EXPECT_EQ(gaps.size(), 2);

    size_t l = 3;
    size_t r = 5;

    for (const auto& i : gaps)
    {
        EXPECT_EQ(i.Left(),  l);
        EXPECT_EQ(i.Right(), r);

        l += 4;
        r += 4;
    }
}

TEST(IntervalTest, Gaps2)
{
    IntervalTree tree;

    tree.Insert(Interval(3, 9));

    IntervalTree gaps = tree.Gaps(Interval(5, 11));

    EXPECT_EQ(gaps.size(), 1);

    for (auto& i : gaps)
    {
        EXPECT_EQ(i.Left(),  9);
        EXPECT_EQ(i.Right(), 11);
    }

    gaps = tree.Gaps(Interval(1, 11));

    EXPECT_EQ(gaps.size(), 2);

    size_t l = 1;
    size_t r = 3;

    for (const auto& i : gaps)
    {
        EXPECT_EQ(i.Left(),  l);
        EXPECT_EQ(i.Right(), r);

        l += 8;
        r += 8;
    }

    gaps = tree.Gaps(Interval(11, 15));

    EXPECT_EQ(gaps.size(), 1);

    for (const auto& i : gaps)
    {
        EXPECT_EQ(i.Left(),  11);
        EXPECT_EQ(i.Right(), 15);
    }
}

TEST(IntervalTest, Gaps3)
{
    IntervalTree tree;

    tree.Insert(Interval(3, 5));
    tree.Insert(Interval(7, 9));

    IntervalTree gaps = tree.Gaps(Interval(4, 9));

    EXPECT_EQ(gaps.size(), 1);

    for (const auto& i : gaps)
    {
        EXPECT_EQ(i.Left(),  5);
        EXPECT_EQ(i.Right(), 7);
    }
}

TEST(IntervalTest, ZMW25300)
{
    IntervalTree tree;

    tree.Insert(Interval(252, 295));
    tree.Insert(Interval(293, 338));

    for (const auto& i : tree)
    {
        EXPECT_EQ(i.Left(),  252);
        EXPECT_EQ(i.Right(), 338);
    }
}
