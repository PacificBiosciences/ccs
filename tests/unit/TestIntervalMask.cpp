// Copyright (c) 2017-2017, Pacific Biosciences of California, Inc.
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
#include <stdexcept>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <pacbio/data/Read.h>

#include <pacbio/consensus/IntervalMask.h>
#include <pacbio/consensus/Mutation.h>

using namespace PacBio::Consensus;

TEST(IntervalMaskTest, Left)
{
    IntervalMask mask;

    mask.Insert({2, 4});

    mask.Mutate({Mutation::Deletion(1, 1)});

    EXPECT_EQ(mask.size(), 1);
    for (const auto& i : mask) {
        EXPECT_EQ(i.Left(), 1);
        EXPECT_EQ(i.Right(), 3);
    }

    mask.Mutate({Mutation::Insertion(1, 'A')});

    EXPECT_EQ(mask.size(), 1);
    for (const auto& i : mask) {
        EXPECT_EQ(i.Left(), 2);
        EXPECT_EQ(i.Right(), 4);
    }

    mask.Mutate({Mutation::Substitution(1, 'A')});

    EXPECT_EQ(mask.size(), 1);
    for (const auto& i : mask) {
        EXPECT_EQ(i.Left(), 2);
        EXPECT_EQ(i.Right(), 4);
    }
}

TEST(IntervalMaskTest, Inside)
{
    IntervalMask mask;

    mask.Insert({2, 5});

    mask.Mutate({Mutation::Deletion(3, 1)});

    EXPECT_EQ(mask.size(), 1);
    for (const auto& i : mask) {
        EXPECT_EQ(i.Left(), 2);
        EXPECT_EQ(i.Right(), 4);
    }

    mask.Mutate({Mutation::Insertion(3, 'A')});

    EXPECT_EQ(mask.size(), 1);
    for (const auto& i : mask) {
        EXPECT_EQ(i.Left(), 2);
        EXPECT_EQ(i.Right(), 5);
    }

    mask.Mutate({Mutation::Substitution(3, 'A')});

    EXPECT_EQ(mask.size(), 1);
    for (const auto& i : mask) {
        EXPECT_EQ(i.Left(), 2);
        EXPECT_EQ(i.Right(), 5);
    }
}

TEST(IntervalMaskTest, Right)
{
    IntervalMask mask;

    mask.Insert({2, 5});

    mask.Mutate({Mutation::Deletion(5, 1)});

    EXPECT_EQ(mask.size(), 1);
    for (const auto& i : mask) {
        EXPECT_EQ(i.Left(), 2);
        EXPECT_EQ(i.Right(), 5);
    }

    mask.Mutate({Mutation::Insertion(5, 'A')});

    EXPECT_EQ(mask.size(), 1);
    for (const auto& i : mask) {
        EXPECT_EQ(i.Left(), 2);
        EXPECT_EQ(i.Right(), 5);
    }

    mask.Mutate({Mutation::Substitution(5, 'A')});

    EXPECT_EQ(mask.size(), 1);
    for (const auto& i : mask) {
        EXPECT_EQ(i.Left(), 2);
        EXPECT_EQ(i.Right(), 5);
    }
}

TEST(IntervalMaskTest, Delete)
{
    IntervalMask mask;

    mask.Insert({1, 2});

    EXPECT_EQ(mask.size(), 1);

    mask.Mutate({Mutation::Deletion(1, 1)});

    EXPECT_EQ(mask.size(), 0);
}

TEST(IntervalMaskTest, Complex)
{
    using PacBio::Data::Interval;

    IntervalMask mask;

    mask.Insert({3, 5});
    mask.Insert({5, 6});  // overlaps, now 3--6
    mask.Insert({9, 12});

    EXPECT_EQ(mask.size(), 2);

    mask.Mutate({Mutation::Insertion(3, 'A'),     // {4, 7}, {10, 13}
                 Mutation::Deletion(4, 1),        // {4, 6}, {9, 12}
                 Mutation::Deletion(8, 1),        // {4, 6}, {8, 11}
                 Mutation::Insertion(10, 'A'),    // {4, 6}, {8, 12}
                 Mutation::Insertion(12, 'A')});  // {4, 6}, {8, 12}

    const std::vector<Interval> ivals(mask.begin(), mask.end());
    const std::vector<Interval> truth{{4, 6}, {8, 12}};

    EXPECT_THAT(ivals, ::testing::ContainerEq(truth));

    mask.Mutate({Mutation::Deletion(4, 1), Mutation::Deletion(5, 1)});

    EXPECT_EQ(mask.size(), 1);
    for (const auto& i : mask) {
        EXPECT_EQ(i.Left(), 6);
        EXPECT_EQ(i.Right(), 10);
    }
}

TEST(IntervalMaskTest, ContainsMutations)
{
    IntervalMask mask;

    mask.Insert({3, 6});

    EXPECT_FALSE(mask.Contains(Mutation::Insertion(3, 'A')));
    EXPECT_TRUE(mask.Contains(Mutation::Insertion(4, 'A')));
    EXPECT_TRUE(mask.Contains(Mutation::Insertion(5, 'A')));
    EXPECT_FALSE(mask.Contains(Mutation::Insertion(6, 'A')));

    EXPECT_FALSE(mask.Contains(Mutation::Deletion(2, 1)));
    EXPECT_TRUE(mask.Contains(Mutation::Deletion(3, 1)));
    EXPECT_TRUE(mask.Contains(Mutation::Deletion(4, 1)));
    EXPECT_TRUE(mask.Contains(Mutation::Deletion(5, 1)));
    EXPECT_FALSE(mask.Contains(Mutation::Deletion(6, 1)));

    EXPECT_FALSE(mask.Contains(Mutation::Substitution(2, 'A')));
    EXPECT_TRUE(mask.Contains(Mutation::Substitution(3, 'A')));
    EXPECT_TRUE(mask.Contains(Mutation::Substitution(4, 'A')));
    EXPECT_TRUE(mask.Contains(Mutation::Substitution(5, 'A')));
    EXPECT_FALSE(mask.Contains(Mutation::Substitution(6, 'A')));
}
