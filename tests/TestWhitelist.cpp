// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
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

#include <pacbio/ccs/Whitelist.h>

using namespace PacBio::CCS;

TEST(WhitelistTest, AllTest)
{
    Whitelist wl1("all"), wl2("*:*");
    EXPECT_TRUE(wl1.Contains("movieName", 34) && wl2.Contains("movieName", 42));
}

TEST(WhitelistTest, CrazyTests)
{
    EXPECT_THROW(Whitelist("1-3;movieName:*"), std::invalid_argument);
    EXPECT_THROW(Whitelist("movieName:*;1-3"), std::invalid_argument);
    EXPECT_THROW(Whitelist("all;1-3"), std::invalid_argument);
    EXPECT_THROW(Whitelist("1-3;all"), std::invalid_argument);
    EXPECT_THROW(Whitelist("movieName:1-3;movieName:4-5"), std::invalid_argument);
}

TEST(WhitelistTest, SingleRange)
{
    Whitelist wl("1-3");

    EXPECT_TRUE(wl.Contains("", 1) && wl.Contains("", 2) && wl.Contains("", 3));
    EXPECT_FALSE(wl.Contains("", 0) || wl.Contains("", 4));
}

TEST(WhitelistTest, TwoMovieRanges)
{
    Whitelist wl("movie1:*;movie2:1-3");

    EXPECT_TRUE(wl.Contains("movie1", 42));
    EXPECT_TRUE(wl.Contains("movie2", 3));
    EXPECT_FALSE(wl.Contains("movie2", 4));
    EXPECT_FALSE(wl.Contains("movie3", 1));
}
