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
    Whitelist wls[] = {Whitelist("1-3"), Whitelist("*:1-3")};

    for (const auto& wl : wls) {
        EXPECT_TRUE(wl.Contains("", 1) && wl.Contains("", 2) && wl.Contains("", 3));
        EXPECT_FALSE(wl.Contains("", 0) || wl.Contains("", 4));
    }
}

TEST(WhitelistTest, TwoMovieRanges)
{
    Whitelist wl("movie1:*;movie2:1-3");

    EXPECT_TRUE(wl.Contains("movie1", 42));
    EXPECT_TRUE(wl.Contains("movie2", 3));
    EXPECT_FALSE(wl.Contains("movie2", 4));
    EXPECT_FALSE(wl.Contains("movie3", 1));
}
