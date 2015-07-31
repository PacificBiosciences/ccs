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

// Author: David Alexander

#include <iostream>

#include <boost/assign.hpp>
#include <boost/assign/std/vector.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <ConsensusCore/Poa/PoaConsensus.hpp>

#include <pacbio/ccs/SparsePoa.h>

#include "TestData.h"
#include "TestUtils.h"

using std::vector;
using std::string;

using namespace boost::assign;

using namespace PacBio::CCS;
using ::testing::Ge;

TEST(SparsePoaTest, TestLocalStaggered)
{
    // Adapted from Pat's C# test
    vector<std::string> reads;

    //        0123456789012345678901234567890
    reads += "TTTACAGGATAGTGCCGCCAATCTTCCAGT",
    //               0123456789012345678901234567890123456789012345
                    "GATACCCCGTGCCGCCAATCTTCCAGTATATACAGCACGAGTAGC",
    //                012345678901234567890123456789012345678901234567890123456789012345678
                     "ATAGTGCCGCCAATCTTCCAGTATATACAGCACGGAGTAGCATCACGTACGTACGTCTACACGTAATT",
    //                                                                    0123456789012345678901234567890123456
                                                                         "ACGTCTACACGTAATTTTGGAGAGCCCTCTCTCACG",
    //                                                                          01234567890123456789012345678901
                                                                               "ACACGTAATTTTGGAGAGCCCTCTCTTCACG",
    //             01234567890123456789012345678901234567890123456789012345
                  "AGGATAGTGCCGCCAATCTTCCAGTAATATACAGCACGGAGTAGCATCACGTACG",
    //                01234567890123456789012345678901234567890123456789012345678901234
                     "ATAGTGCCGCCAATCTTCCAGTATATACAGCACGGAGTAGCATCACGTACGTACGTCTACACGT";
    // -----------------------------------------------------------------------------------
    //                012345678901234567890123456789012345678901234567890123456789012345678
    //               "ATAGTGCCGCCAATCTTCCAGTATATACAGCACGGAGTAGCATCACGTACGTACGTCTACACGTAATT",

    SparsePoa sp;
    for (auto& read : reads)
    {
        SparsePoa::ReadKey id = sp.OrientAndAddRead(read);
        EXPECT_THAT(id, Ge(0));
    }

    vector<PoaAlignmentSummary> summaries;
    string consensusSeq = sp.FindConsensus(4, &summaries)->Sequence;

    EXPECT_EQ("ATAGTGCCGCCAATCTTCCAGTATATACAGCACGGAGTAGCATCACGTACGTACGTCTACACGTAATT", consensusSeq);

    EXPECT_FALSE(summaries[0].ReverseComplementedRead);
    EXPECT_EQ(Interval( 8, 30), summaries[0].ExtentOnRead);
    EXPECT_EQ(Interval( 0, 22), summaries[0].ExtentOnConsensus);

    EXPECT_FALSE(summaries[1].ReverseComplementedRead);
    EXPECT_EQ(Interval( 8, 45), summaries[1].ExtentOnRead);
    EXPECT_EQ(Interval( 3, 41), summaries[1].ExtentOnConsensus);

    EXPECT_FALSE(summaries[2].ReverseComplementedRead);
    EXPECT_EQ(Interval( 0, 68), summaries[2].ExtentOnRead);
    EXPECT_EQ(Interval( 0, 68), summaries[2].ExtentOnConsensus);

    EXPECT_FALSE(summaries[3].ReverseComplementedRead);
    EXPECT_EQ(Interval( 0, 16), summaries[3].ExtentOnRead);
    EXPECT_EQ(Interval(52, 68), summaries[3].ExtentOnConsensus);

    EXPECT_FALSE(summaries[4].ReverseComplementedRead);
    EXPECT_EQ(Interval( 0, 10), summaries[4].ExtentOnRead);
    EXPECT_EQ(Interval(58, 68), summaries[4].ExtentOnConsensus);

    EXPECT_FALSE(summaries[5].ReverseComplementedRead);
    EXPECT_EQ(Interval( 3, 55), summaries[5].ExtentOnRead);
    EXPECT_EQ(Interval( 0, 51), summaries[5].ExtentOnConsensus);

    EXPECT_FALSE(summaries[6].ReverseComplementedRead);
    EXPECT_EQ(Interval(0, 64),  summaries[6].ExtentOnRead);
    EXPECT_EQ(Interval(0, 64),  summaries[6].ExtentOnConsensus);
}


TEST(SparsePoaTest, TestOrientation)
{
    vector<std::string> reads = { "AAAGATTACAGGG",
                                  "CCCTGTAATCTTT",
                                  "AAAGATTACAGGG" };
    SparsePoa sp;
    for (auto& read : reads)
    {
        SparsePoa::ReadKey id = sp.OrientAndAddRead(read);
        EXPECT_THAT(id, Ge(0));
    }


    vector<PoaAlignmentSummary> summaries;
    string consensusSeq = sp.FindConsensus(2, &summaries)->Sequence;

    EXPECT_EQ("AAAGATTACAGGG", consensusSeq);

    EXPECT_FALSE(summaries[0].ReverseComplementedRead);
    EXPECT_TRUE (summaries[1].ReverseComplementedRead);
    EXPECT_FALSE(summaries[2].ReverseComplementedRead);
}


TEST(SparsePoa, TestZmw6251)
{
    using std::cout;
    using std::endl;

    std::string fastaFname = tests::DataDir + "/m140905_042212_sidney_c100564852550000001823085912221377_s1_X0.fasta";
    vector<string> ids, seqs;
    LoadFastaSequences(fastaFname, ids, seqs);

    SparsePoa sp;
    for (const auto& seq : seqs)
    {
        SparsePoa::ReadKey id = sp.OrientAndAddRead(seq);
        EXPECT_THAT(id, Ge(0));
    }

    vector<PoaAlignmentSummary> summaries;
    auto pc = sp.FindConsensus(8, &summaries);
    ASSERT_EQ(10, pc->Graph.NumReads());

    string consensusSeq = pc->Sequence;
    //pc->WriteGraphVizFile("/tmp/zmw6251.dot");

    // What it looks like:
    //
    // r0:     >>>>>>>>>>>
    // r1: <<<<<<<<<<<<<<<
    // r2: >>>>>>>>>>>>>>>
    // ..
    // r8: >>>>>>>>>>>>>>>
    // r9:           <<<<<
    for (int i=0; i<=9; i++)
    {
        if (i % 2 == 0)
            EXPECT_FALSE(summaries[0].ReverseComplementedRead);
        else
            EXPECT_TRUE (summaries[1].ReverseComplementedRead);
    }

    // css ~ 600 bases; check that things hit roughly as expected
    EXPECT_TRUE(summaries[0].ExtentOnConsensus.Covers(Interval(300, 595)));
    for (int i=1; i <= 8; i++)
        EXPECT_TRUE(summaries[1].ExtentOnConsensus.Covers(Interval(5, 595)));
    EXPECT_TRUE(summaries[0].ExtentOnConsensus.Covers(Interval(300, 595)));
}
