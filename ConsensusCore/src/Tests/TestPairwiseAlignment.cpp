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
#include <boost/shared_ptr.hpp>

#include <ConsensusCore/Align/AffineAlignment.hpp>
#include <ConsensusCore/Align/LinearAlignment.hpp>
#include <ConsensusCore/Align/PairwiseAlignment.hpp>


using namespace ConsensusCore;  // NOLINT
using ::testing::ElementsAreArray;

TEST(PairwiseAlignmentTests, RepresentationTests)
{
    PairwiseAlignment a("GATC", "GA-C");
    EXPECT_EQ("GATC", a.Target());
    EXPECT_EQ("GA-C", a.Query());
    EXPECT_EQ(4, a.Length());
    EXPECT_EQ(3, a.Matches());
    EXPECT_EQ(1, a.Deletions());
    EXPECT_EQ(0, a.Mismatches());
    EXPECT_EQ(0, a.Insertions());
    EXPECT_FLOAT_EQ(0.75, a.Accuracy());
    EXPECT_EQ("MMDM", a.Transcript());

    PairwiseAlignment a2("GATTA-CA",
                         "CA-TAACA");
    EXPECT_EQ("RMDMMIMM", a2.Transcript());
    EXPECT_FLOAT_EQ(5./8, a2.Accuracy());
    EXPECT_EQ(1, a2.Mismatches());
    EXPECT_EQ(1, a2.Deletions());
    EXPECT_EQ(1, a2.Insertions());
    EXPECT_EQ(5, a2.Matches());
}


TEST(PairwiseAlignmentTests, GlobalAlignmentTests)
{
    PairwiseAlignment* a = Align("GATT", "GATT");
    EXPECT_FLOAT_EQ(1.0, a->Accuracy());
    EXPECT_EQ("GATT", a->Target());
    EXPECT_EQ("GATT", a->Query());
    EXPECT_EQ("MMMM", a->Transcript());
    delete a;

    a = Align("GATT", "GAT");
    EXPECT_FLOAT_EQ(0.75, a->Accuracy());
    EXPECT_EQ("GATT", a->Target());
    EXPECT_EQ("GA-T", a->Query());
    EXPECT_EQ("MMDM", a->Transcript());
    delete a;

    a = Align("GATTACA", "TT");
    EXPECT_EQ("GATTACA", a->Target());
    EXPECT_EQ("--TT---", a->Query());
    EXPECT_FLOAT_EQ(2./7, a->Accuracy());
    delete a;
}


TEST(PairwiseAlignmentTests, TargetPositionsInQueryTest)
{
    // MMM -> 0123
    {
        int expected[] = { 0, 1, 2, 3 };
        ASSERT_THAT(TargetToQueryPositions("MMM"), ElementsAreArray(expected));
    }

    // DMM -> 0012, MDM -> 0112, MMD -> 0122,
    {
        int expected1[] = { 0, 0, 1, 2 };
        int expected2[] = { 0, 1, 1, 2 };
        int expected3[] = { 0, 1, 2, 2 };
        ASSERT_THAT(TargetToQueryPositions("DMM"), ElementsAreArray(expected1));
        ASSERT_THAT(TargetToQueryPositions("MDM"), ElementsAreArray(expected2));
        ASSERT_THAT(TargetToQueryPositions("MMD"), ElementsAreArray(expected3));
    }

    // IMM -> 123, MIM -> 023, MMI -> 013,
    {
        int expected1[] = { 1, 2, 3 };
        int expected2[] = { 0, 2, 3 };
        int expected3[] = { 0, 1, 3 };
        ASSERT_THAT(TargetToQueryPositions("IMM"), ElementsAreArray(expected1));
        ASSERT_THAT(TargetToQueryPositions("MIM"), ElementsAreArray(expected2));
        ASSERT_THAT(TargetToQueryPositions("MMI"), ElementsAreArray(expected3));
    }

    // MRM, MDIM -> 0123
    // MIDM -> 0223
    {
        int expected1[] = { 0, 1, 2, 3 };
        int expected2[] = { 0, 2, 2, 3 };
        ASSERT_THAT(TargetToQueryPositions("MRM"),  ElementsAreArray(expected1));
        ASSERT_THAT(TargetToQueryPositions("MDIM"), ElementsAreArray(expected1));
        ASSERT_THAT(TargetToQueryPositions("MIDM"), ElementsAreArray(expected2));
    }
}


// ------------------ AffineAlignment tests ---------------------

TEST(AffineAlignmentTests, BasicTests)
{
    PairwiseAlignment* a = AlignAffine("ATT", "ATT");
    EXPECT_EQ("ATT", a->Target());
    EXPECT_EQ("ATT", a->Query());
    delete a;

    a = AlignAffine("AT", "ATT");
    EXPECT_EQ("A-T", a->Target());
    EXPECT_EQ("ATT", a->Query());
    delete a;

    a = AlignAffine("GA", "GAT");
    EXPECT_EQ("GA-", a->Target());
    EXPECT_EQ("GAT", a->Query());
    delete a;

    a = AlignAffine("GAT", "GA");
    EXPECT_EQ("GAT", a->Target());
    EXPECT_EQ("GA-", a->Query());
    delete a;

    a = AlignAffine("GA", "TGA");
    EXPECT_EQ("-GA", a->Target());
    EXPECT_EQ("TGA", a->Query());
    delete a;

    a = AlignAffine("TGA", "GA");
    EXPECT_EQ("TGA", a->Target());
    EXPECT_EQ("-GA", a->Query());
    delete a;

    a = AlignAffine("GATTACA", "GATTTACA");
    EXPECT_EQ("GA-TTACA", a->Target());
    EXPECT_EQ("GATTTACA", a->Query());
    delete a;
}


TEST(AffineAlignmentTests, LargeGapTest)
{
    // Test a real-world large insertion, found in an E. Coli
    // experiment
    const char* target = \
        "AACGATTTTATGATGGCATGTGACATGTATTTCCGTTGGGGGCATTTTAATAAGTGAGGA"
        "AGTGATAGGAAGTGACCAGATAATACATATATGTTCTGTACTCTCTTGCGCATTTTGATT"
        "GTTGACTGAGTAACCAGACAGTTGATGTGCACGATTTCCCCTCGCCCTAACAGACGTGGG"
        "CGGGGGCACCCGGTGTGTGCAATATACGAAAAAAAAGCCCGTACTTTCGTACGAGCTCTT"
        "CTTTAAATATGGCGGTGAGGGGGGGATTGACTCGCTTCGCTCGCCCTGCGGGCAGCCCGC"
        "TCACTGCGTTCACGGTCTGTCCAACTGGCTGTCGCCAGTTGTCGACCCCCGGTCGGGGCT"
        "TCTCATCCCCCCGGTGTGTGCAATACACGAAAAAAAAGCCCGTACTTTCGTACGAGCTCT"
        "TCTTTAAATATGGCGGTGAGGGGGGGATTGACTCGCTTCGCTCGCCCTGCGGGCAGCCCG"
        "CTCACTGCGTTCACGGTCTGTCCAACTGGCTGTCGCCAGTTGTCGAACCCCGGTCGGGGC"
        "TTCTCATCCCCCCGGTGTGTGCAATATACGAAAAAAAAGCCCGTACTTTCGTACGAGCTC"
        "TTCTTTAAATATGGCGGTGAGGGGGGGATTGACTCGCTTCGCTCGCCCTGCGGGCAGCCC"
        "GCTCACTGCGTTCACGGTCTGTCCAACTGGCTGTCGCCAGTTGTCGAACCCCGGTCGGGG"
        "CTTCTCATCCCCCCGGTGTGTGCAATATACGAAAAAAAAGCCCGTACTTTCGTACGAGCT"
        "CTTCTTTAAATATGGCGGTGAGGGGGGGATTCGAACCCCCGATACGTTGCCGTATACACA"
        "CTTTCCAGGCGTGCTCCTTCAGCCACTCGGACACCTCACCAAATTGTCGTTCCTGTCTTG"
        "CTGGAACGGGCGCTAATTTAGGGAAATCATGACCTGAGGTCAACAAACTTTTTGAAAAAA"
        "TCGCGCGTTTATTCAAACTTCAATCAATGTGTGGTTTTAATAAGCGAAAT";

    const char* query = \
        "AACGATTTTATGATGGCATGTGACATGTATTTCCGTTGGGGGCATTTTAATAAGTGAGGA"
        "AGTGATAGGAAGTGACCAGATAATACATATATGTTCTGTACTCTCTTGCGCATTTTGATT"
        "GTTGACTGAGTAACCAGACAGTTGATGTGCACGATTTCCCCTCGCCCTAACAGACGTGGG"
        "CGGGGGCACCCGGTGTGTGCAATATACGAAAAAAAAGCCCGTACTTTCGTACGAGCTCTT"
        "CTTTAAATATGGCGGTGAGGGGGGGATTGACTCGCTTCGCTCGCCCTGCGGGCAGCCCGC"
        "TCACTGCGTTCACGGTCTGTCCAACTGGCTGTCGCCAGTTGTCGACCCCCGGTCGGGGCT"
        "TCTCATCCCCCCGGTGTGTGCAATACACGAAAAAAAAGCCCGTACTTTCGTACGAGCTCT"
        "TCTTTAAATATGGCGGTGAGGGGGGGATTGACTCGCTTCGCTCGCCCTGCGGGCAGCCCG"
        "CTCACTGCGTTCACGGTCTGTCCAACTGGCTGTCGCCAGTTGTCGAACCCCGGTCGGGGC"
        "TTCTCATCCCCCCGGTGTGTGCAATATACGAAAAAAAAGCCCGTACTTTCGTACGAGCTC"
        "TTCTTTAAATATGGCGGTGAGGGGGGGATTGACTCGCTTCGCTCGCCCTGCGGGCAGCCC"
        "GCTCACTGCGTTCACGGTCTGTCCAACTGGCTGTCGCCAGTTGTCGAACCCCGGTCGGGG"
        "CTTCTCATCCCCCCGGTGTGTGCAATATACGAAAAAAAAGCCCGTACTTTCGTACGAGCT"
        "CTTCTTTAAATATGGCGGTGAGGGGGGGATTGACTCGCTTCGCTCGCCCTGCGGGCAGCC"
        "CGCTCACTGCGTTCACGGTCTGTCCAACTGGCTGTCGCCAGTTGTCGAACCCCGGTCGGG"
        "GCTTCTCATCCCCCCGGTGTGTGCAATATACGAAAAAAAAGCCCGTACTTTCGTACGAGC"
        "TCTTCTTTAAATATGGCGGTGAGGGGGGGATTCGAACCCCCGATACGTTGCCGTATACAC"
        "ACTTTCCAGGCGTGCTCCTTCAGCCACTCGGACACCTCACCAAATTGTCGTTCCTGTCTT"
        "GCTGGAACGGGCGCTAATTTAGGGAAATCATGACCTGAGGTCAACAAACTTTTTGAAAAA"
        "ATCGCGCGTTTATTCAAACTTCAATCAATGTGTGGTTTTAATAAGCGAAAT";

    const char* expectedAlignedTarget = \
        "AACGATTTTATGATGGCATGTGACATGTATTTCCGTTGGGGGCATTTTAATAAGTGAGGA"
        "AGTGATAGGAAGTGACCAGATAATACATATATGTTCTGTACTCTCTTGCGCATTTTGATT"
        "GTTGACTGAGTAACCAGACAGTTGATGTGCACGATTTCCCCTCGCCCTAACAGACGTGGG"
        "CGGGGGCACCCGGTGTGTGCAATATACGAAAAAAAAGCCCGTACTTTCGTACGAGCTCTT"
        "CTTTAAATATGGCGGTGAGGGGGGGATTGACTCGCTTCGCTCGCCCTGCGGGCAGCCCGC"
        "TCACTGCGTTCACGGTCTGTCCAACTGGCTGTCGCCAGTTGTCGACCCCCGGTCGGGGCT"
        "TCTCATCCCCCCGGTGTGTGCAATAC----------------------------------"
        "------------------------------------------------------------"
        "------------------------------------------------------------"
        "---------------------------ACGAAAAAAAAGCCCGTACTTTCGTACGAGCTC"
        "TTCTTTAAATATGGCGGTGAGGGGGGGATTGACTCGCTTCGCTCGCCCTGCGGGCAGCCC"
        "GCTCACTGCGTTCACGGTCTGTCCAACTGGCTGTCGCCAGTTGTCGAACCCCGGTCGGGG"
        "CTTCTCATCCCCCCGGTGTGTGCAATATACGAAAAAAAAGCCCGTACTTTCGTACGAGCT"
        "CTTCTTTAAATATGGCGGTGAGGGGGGGATTGACTCGCTTCGCTCGCCCTGCGGGCAGCC"
        "CGCTCACTGCGTTCACGGTCTGTCCAACTGGCTGTCGCCAGTTGTCGAACCCCGGTCGGG"
        "GCTTCTCATCCCCCCGGTGTGTGCAATATACGAAAAAAAAGCCCGTACTTTCGTACGAGC"
        "TCTTCTTTAAATATGGCGGTGAGGGGGGGATTCGAACCCCCGATACGTTGCCGTATACAC"
        "ACTTTCCAGGCGTGCTCCTTCAGCCACTCGGACACCTCACCAAATTGTCGTTCCTGTCTT"
        "GCTGGAACGGGCGCTAATTTAGGGAAATCATGACCTGAGGTCAACAAACTTTTTGAAAAA"
        "ATCGCGCGTTTATTCAAACTTCAATCAATGTGTGGTTTTAATAAGCGAAAT";

    PairwiseAlignment* a = AlignAffine(target, query);
    ASSERT_EQ(expectedAlignedTarget, a->Target());
    delete a;
}



// ------------------ IUPAC-aware alignment tests ---------------------

TEST(IupacAlignmentTests, BasicTest)
{
    PairwiseAlignment* a;
    a = AlignAffineIupac("GATTTT", "GMTTT");
    ASSERT_EQ("GATTTT", a->Target());
    ASSERT_EQ("GM-TTT", a->Query());
    delete a;

    a = AlignAffineIupac("TTTTAG", "TTTMG");
    ASSERT_EQ("TTTTAG", a->Target());
    ASSERT_EQ("-TTTMG", a->Query());
    delete a;
}


// ---------------- Linear-space alignment tests -----------------------

TEST(LinearAlignmentTests, BasicTest)
{
    AlignParams params(2, -1, -2, -2);
    AlignConfig config(params, GLOBAL);

    int score, peerScore;
    PairwiseAlignment *a, *peerAlignment;

    a = AlignLinear("GATTACA", "GATTACA", &score);
    EXPECT_EQ("GATTACA", a->Target());
    EXPECT_EQ("GATTACA", a->Query());
    EXPECT_EQ("MMMMMMM", a->Transcript());
    EXPECT_EQ(14, score);
    delete a;

    a = AlignLinear("TATGC", "AGTACGCA", &score);
    EXPECT_EQ("--TATGC-", a->Target());
    EXPECT_EQ("AGTACGCA", a->Query());
    EXPECT_EQ("IIMMRMMI", a->Transcript());
    EXPECT_EQ(1, score);
    delete a;

    a = AlignLinear("AGTACGCA", "TATGC", &score);
    EXPECT_EQ("AGTACGCA", a->Target());
    EXPECT_EQ("--TATGC-", a->Query());
    EXPECT_EQ("DDMMRMMD", a->Transcript());
    EXPECT_EQ(1, score);
    delete a;


    a = AlignLinear("GATT", "GATT");
    EXPECT_FLOAT_EQ(1.0, a->Accuracy());
    EXPECT_EQ("GATT", a->Target());
    EXPECT_EQ("GATT", a->Query());
    EXPECT_EQ("MMMM", a->Transcript());
    delete a;

    a = AlignLinear("GATT", "GAT");
    EXPECT_FLOAT_EQ(0.75, a->Accuracy());
    EXPECT_EQ("GATT", a->Target());
    EXPECT_EQ("GA-T", a->Query());
    EXPECT_EQ("MMDM", a->Transcript());
    delete a;

    a = AlignLinear("GATTACA", "TT");
    EXPECT_EQ("GATTACA", a->Target());
    EXPECT_EQ("--TT---", a->Query());
    EXPECT_FLOAT_EQ(2./7, a->Accuracy());
    delete a;

    const char* ref = "GTATTTTAAATAAAAACATTAAGTTATGACGAAGAAGAACGGAAACGCCTTAAACCGGAAAATTTTCATAAATAGCGAAAACCCGCGAGGTCGCCGCCC";
    const char* read = "GTATTTTAAATAAAAAAACATTATAGTTTAATGAACGAGAATGAACGGTAATACGCCTTTAAAGCCTGAAATATTTTTCCATAAATGTAATTTCTGTATATAATCTCCGCGAGTGTCTGCCGCCC";

    a = AlignLinear(ref, read, &score);
    peerAlignment = Align(ref, read, &peerScore, config);
    EXPECT_EQ(score, peerScore);
}


#if 0
TEST(LinearAlignmentTests, SemiglobalTests)
{
    AlignParams params(2, -1, -2, -2);
    int score, peerScore;
    PairwiseAlignment *a, *peerAlignment;

    a = AlignLinear("AGTCGATACACCCC", "GATTACA");
    EXPECT_EQ("AGTCGA-TACACCCC", a->Target());
    EXPECT_EQ("----GATTACA----", a->Query());

    // we got:
    // -AGTCGATACACCCC
    // GA-T---TACA----

}
#endif
