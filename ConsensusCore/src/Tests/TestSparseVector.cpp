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

#include <gtest/gtest.h>

#include <boost/format.hpp>
#include <iostream>
#include <string>
#include <vector>

#include <ConsensusCore/Matrix/SparseVector.hpp>

using namespace ConsensusCore; // NOLINT

TEST(SparseVectorTest, BasicTest)
{
    SparseVectorF sv(100, 10, 20);
    EXPECT_LE(10, sv.AllocatedEntries());

    for (int i = 0; i < 100; i++)
    {
        EXPECT_EQ(sv(i), -FLT_MAX);
    }

    for (int i = 10; i < 20; i++)
    {
        sv.Set(i, i);
    }
    for (int i = 0; i < 100; i++)
    {
        if (i >= 10 && i < 20) EXPECT_EQ(i, sv(i));
        else EXPECT_EQ(-FLT_MAX, sv(i)); // NOLINT
    }

    sv.Set(50, 50);
    EXPECT_LE(40, sv.AllocatedEntries());
    for (int i = 0; i < 100; i++)
    {
        if (i >= 10 && i < 20) EXPECT_EQ(i, sv(i));
        else if (i == 50) EXPECT_EQ(i, sv(i));
        else EXPECT_EQ(-FLT_MAX, sv(i)); // NOLINT
    }
}


TEST(SparseVectorTest, BasicTest2)
{
    SparseVectorF sv(100, 50, 60);

    sv.Set(5, 5);
    for (int i = 0; i < 100; i++)
    {
        if (i == 5) EXPECT_EQ(i, sv(i));
        else EXPECT_EQ(-FLT_MAX, sv(i)); // NOLINT
    }
}



TEST(SparseVector, CopyTest)
{
    SparseVectorF sv(10, 3, 7);
    sv.Set(4, 5);

    SparseVectorF svCopy(sv);
    ASSERT_EQ(5, svCopy(4));

    for (int i = 0; i < 10; i++)
    {
        ASSERT_EQ(sv(i), svCopy(i));
    }
}
