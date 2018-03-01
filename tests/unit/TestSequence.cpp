// Author: Lance Hepler

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <string>

#include <pacbio/data/Sequence.h>

using std::string;

using namespace PacBio::Data;  // NOLINT

namespace SequenceTests {

TEST(SequenceTest, ReverseComplement) { EXPECT_EQ("TACGAT", ReverseComplement("ATCGTA")); }

TEST(SequenceTest, ReverseComplementAmbiguous)
{
    EXPECT_EQ("N-YRWSKMDVHB-N", ReverseComplement("N-VDBHKMSWYR-N"));
}

TEST(SequenceTest, ReverseComplementInvalid)
{
    // test space of invalid letters exhaustively
    EXPECT_ANY_THROW(ReverseComplement("E"));
    EXPECT_ANY_THROW(ReverseComplement("F"));
    EXPECT_ANY_THROW(ReverseComplement("I"));
    EXPECT_ANY_THROW(ReverseComplement("J"));
    EXPECT_ANY_THROW(ReverseComplement("L"));
    EXPECT_ANY_THROW(ReverseComplement("O"));
    EXPECT_ANY_THROW(ReverseComplement("P"));
    EXPECT_ANY_THROW(ReverseComplement("Q"));
    EXPECT_ANY_THROW(ReverseComplement("U"));
    EXPECT_ANY_THROW(ReverseComplement("X"));
    EXPECT_ANY_THROW(ReverseComplement("Z"));
}
}
