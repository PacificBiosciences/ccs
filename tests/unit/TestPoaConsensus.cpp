// Author: David Alexander

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/assign.hpp>
#include <boost/assign/std/vector.hpp>
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <pacbio/align/AlignConfig.h>
#include <pacbio/consensus/Mutation.h>
#include <pacbio/denovo/PoaConsensus.h>

using std::string;
using std::vector;
using std::cout;
using std::endl;
using boost::erase_all_copy;

using namespace boost::assign;      // NOLINT
using namespace PacBio::Consensus;  // NOLINT
using namespace PacBio::Poa;        // NOLINT
using namespace PacBio::Align;      // NOLINT

#define MAKE_ALL_PLOTS false

namespace PoaConsensusTests {

static void plotConsensus(const PacBio::Poa::PoaConsensus* pc, string description,
                          bool REALLY_MAKE_THIS_ONE = false)
{
    if (MAKE_ALL_PLOTS || REALLY_MAKE_THIS_ONE) {
        string dotFname = description + ".dot";
        string pngFname = description + ".png";
        string cmd = string("dot -Tpng ") + dotFname + " > " + pngFname;
        pc->Graph.WriteGraphVizFile(description + ".dot",
                                    (PoaGraph::COLOR_NODES | PoaGraph::VERBOSE_NODES), pc);
        // cout << cmd << endl;
        int ret = system(cmd.c_str());
    }
}

}  // namespace PoaConsensusTests

// TEST(PoaGraph, NoReadsTest)
// {
//  // Test that it works with no reads
//  vector<const SequenceFeatures*> reads;
//  const PacBio::Poa::PoaConsensus* pc = PoaConsensus::findConsensus(reads,
//  PoaConfig::GLOBAL_ALIGNMENT);
//  string dot = pc->getGraph()->toGraphViz();
//  cout << dot << endl;
// }

TEST(PoaGraph, SmallBasicTest1)
{
    // Test that it works with a single sequence
    vector<string> reads;
    reads += "GGG";
    const PacBio::Poa::PoaConsensus* pc = PoaConsensus::FindConsensus(reads, AlignMode::GLOBAL);
    string dot = pc->Graph.ToGraphViz();
    string expectedDot =
        "digraph G {"
        "rankdir=\"LR\";"
        "0[shape=Mrecord, label=\"{ ^ | 0 }\"];"
        "1[shape=Mrecord, label=\"{ $ | 0 }\"];"
        "2[shape=Mrecord, label=\"{ G | 1 }\"];"
        "3[shape=Mrecord, label=\"{ G | 1 }\"];"
        "4[shape=Mrecord, label=\"{ G | 1 }\"];"
        "0->2 ;"
        "2->3 ;"
        "3->4 ;"
        "4->1 ;"
        "}";
    PoaConsensusTests::plotConsensus(pc, "small-basic-1");
    EXPECT_EQ(expectedDot, erase_all_copy(dot, "\n"));
    EXPECT_EQ("GGG", pc->Sequence);
    delete pc;
}

TEST(PoaGraph, SmallBasicTest2)
{
    // Test that it works with two identical sequences
    vector<string> reads;
    reads += "GGG", "GGG";
    const PacBio::Poa::PoaConsensus* pc = PoaConsensus::FindConsensus(reads, AlignMode::GLOBAL);
    string dot = pc->Graph.ToGraphViz();
    string expectedDot =
        "digraph G {"
        "rankdir=\"LR\";"
        "0[shape=Mrecord, label=\"{ ^ | 0 }\"];"
        "1[shape=Mrecord, label=\"{ $ | 0 }\"];"
        "2[shape=Mrecord, label=\"{ G | 2 }\"];"
        "3[shape=Mrecord, label=\"{ G | 2 }\"];"
        "4[shape=Mrecord, label=\"{ G | 2 }\"];"
        "0->2 ;"
        "2->3 ;"
        "3->4 ;"
        "4->1 ;"
        "}";
    PoaConsensusTests::plotConsensus(pc, "small-basic-2");
    EXPECT_EQ(expectedDot, erase_all_copy(dot, "\n"));
    EXPECT_EQ("GGG", pc->Sequence);

    delete pc;
}

TEST(PoaGraph, SmallExtraTests)
{
    // Extra at beginning
    {
        vector<string> reads;
        reads += "GGG", "TGGG";
        const PacBio::Poa::PoaConsensus* pc = PoaConsensus::FindConsensus(reads, AlignMode::GLOBAL);
        string dot = pc->Graph.ToGraphViz();
        string expectedDot =
            "digraph G {"
            "rankdir=\"LR\";"
            "0[shape=Mrecord, label=\"{ ^ | 0 }\"];"
            "1[shape=Mrecord, label=\"{ $ | 0 }\"];"
            "2[shape=Mrecord, label=\"{ G | 2 }\"];"
            "3[shape=Mrecord, label=\"{ G | 2 }\"];"
            "4[shape=Mrecord, label=\"{ G | 2 }\"];"
            "5[shape=Mrecord, label=\"{ T | 1 }\"];"
            "0->2 ;"
            "2->3 ;"
            "3->4 ;"
            "4->1 ;"
            "5->2 ;"
            "0->5 ;"
            "}";
        PoaConsensusTests::plotConsensus(pc, "extra-at-beginning");
        EXPECT_EQ(expectedDot, erase_all_copy(dot, "\n"));
        EXPECT_EQ("GGG", pc->Sequence);
        delete pc;
    }

    // Extra in middle
    {
        vector<string> reads;
        reads += "GGG", "GTGG";
        const PacBio::Poa::PoaConsensus* pc = PoaConsensus::FindConsensus(reads, AlignMode::GLOBAL);
        string dot = pc->Graph.ToGraphViz();
        // this would be easier if we could use the C++0x raw
        // strings feature (in g++ 4.5+)
        string expectedDot =
            "digraph G {"
            "rankdir=\"LR\";"
            "0[shape=Mrecord, label=\"{ ^ | 0 }\"];"
            "1[shape=Mrecord, label=\"{ $ | 0 }\"];"
            "2[shape=Mrecord, label=\"{ G | 2 }\"];"
            "3[shape=Mrecord, label=\"{ G | 2 }\"];"
            "4[shape=Mrecord, label=\"{ G | 2 }\"];"
            "5[shape=Mrecord, label=\"{ T | 1 }\"];"
            "0->2 ;"
            "2->3 ;"
            "3->4 ;"
            "4->1 ;"
            "5->3 ;"
            "2->5 ;"
            "}";
        PoaConsensusTests::plotConsensus(pc, "extra-in-middle");
        EXPECT_EQ(expectedDot, erase_all_copy(dot, "\n"));
        EXPECT_EQ("GGG", pc->Sequence);
        delete pc;
    }

    // Extra at end
    {
        vector<string> reads;
        reads += "GGG", "GGGT";
        const PacBio::Poa::PoaConsensus* pc = PoaConsensus::FindConsensus(reads, AlignMode::GLOBAL);
        string dot = pc->Graph.ToGraphViz();
        string expectedDot =
            "digraph G {"
            "rankdir=\"LR\";"
            "0[shape=Mrecord, label=\"{ ^ | 0 }\"];"
            "1[shape=Mrecord, label=\"{ $ | 0 }\"];"
            "2[shape=Mrecord, label=\"{ G | 2 }\"];"
            "3[shape=Mrecord, label=\"{ G | 2 }\"];"
            "4[shape=Mrecord, label=\"{ G | 2 }\"];"
            "5[shape=Mrecord, label=\"{ T | 1 }\"];"
            "0->2 ;"
            "2->3 ;"
            "3->4 ;"
            "4->1 ;"
            "5->1 ;"
            "4->5 ;"
            "}";
        PoaConsensusTests::plotConsensus(pc, "extra-at-end");
        EXPECT_EQ(expectedDot, erase_all_copy(dot, "\n"));
        EXPECT_EQ("GGG", pc->Sequence);
        delete pc;
    }
}

TEST(PoaGraph, SmallMismatchTests)
{
    // Mismatch at beginning
    {
        vector<string> reads;
        reads += "GGG", "TGG";
        const PacBio::Poa::PoaConsensus* pc = PoaConsensus::FindConsensus(reads, AlignMode::GLOBAL);
        string dot = pc->Graph.ToGraphViz();
        string expectedDot =
            "digraph G {"
            "rankdir=\"LR\";"
            "0[shape=Mrecord, label=\"{ ^ | 0 }\"];"
            "1[shape=Mrecord, label=\"{ $ | 0 }\"];"
            "2[shape=Mrecord, label=\"{ G | 1 }\"];"
            "3[shape=Mrecord, label=\"{ G | 2 }\"];"
            "4[shape=Mrecord, label=\"{ G | 2 }\"];"
            "5[shape=Mrecord, label=\"{ T | 1 }\"];"
            "0->2 ;"
            "2->3 ;"
            "3->4 ;"
            "4->1 ;"
            "5->3 ;"
            "0->5 ;"
            "}";
        PoaConsensusTests::plotConsensus(pc, "mismatch-at-beginning");
        EXPECT_EQ(expectedDot, erase_all_copy(dot, "\n"));
        EXPECT_EQ("GG", pc->Sequence);
        delete pc;
    }

    // Mismatch in middle
    {
        vector<string> reads;
        reads += "GGG", "GTG", "GTG";
        const PacBio::Poa::PoaConsensus* pc = PoaConsensus::FindConsensus(reads, AlignMode::GLOBAL);
        string dot = pc->Graph.ToGraphViz();
        string expectedDot =
            "digraph G {"
            "rankdir=\"LR\";"
            "0[shape=Mrecord, label=\"{ ^ | 0 }\"];"
            "1[shape=Mrecord, label=\"{ $ | 0 }\"];"
            "2[shape=Mrecord, label=\"{ G | 3 }\"];"
            "3[shape=Mrecord, label=\"{ G | 1 }\"];"
            "4[shape=Mrecord, label=\"{ G | 3 }\"];"
            "5[shape=Mrecord, label=\"{ T | 2 }\"];"
            "0->2 ;"
            "2->3 ;"
            "3->4 ;"
            "4->1 ;"
            "5->4 ;"
            "2->5 ;"
            "}";
        PoaConsensusTests::plotConsensus(pc, "mismatch-in-middle");
        EXPECT_EQ(expectedDot, erase_all_copy(dot, "\n"));
        EXPECT_EQ("GTG", pc->Sequence);
        delete pc;
    }

    // Mismatch at end
    {
        vector<string> reads;
        reads += "GGG", "GGT";
        const PacBio::Poa::PoaConsensus* pc = PoaConsensus::FindConsensus(reads, AlignMode::GLOBAL);
        string dot = pc->Graph.ToGraphViz();
        string expectedDot =
            "digraph G {"
            "rankdir=\"LR\";"
            "0[shape=Mrecord, label=\"{ ^ | 0 }\"];"
            "1[shape=Mrecord, label=\"{ $ | 0 }\"];"
            "2[shape=Mrecord, label=\"{ G | 2 }\"];"
            "3[shape=Mrecord, label=\"{ G | 2 }\"];"
            "4[shape=Mrecord, label=\"{ G | 1 }\"];"
            "5[shape=Mrecord, label=\"{ T | 1 }\"];"
            "0->2 ;"
            "2->3 ;"
            "3->4 ;"
            "4->1 ;"
            "5->1 ;"
            "3->5 ;"
            "}";
        PoaConsensusTests::plotConsensus(pc, "mismatch-at-end");
        EXPECT_EQ(expectedDot, erase_all_copy(dot, "\n"));
        EXPECT_EQ("GG", pc->Sequence);
        delete pc;
    }
}

TEST(PoaGraph, SmallDeletionTests)
{
    // Deletion at beginning
    {
        vector<string> reads;
        reads += "GAT", "AT";
        const PacBio::Poa::PoaConsensus* pc = PoaConsensus::FindConsensus(reads, AlignMode::GLOBAL);
        string dot = pc->Graph.ToGraphViz();
        string expectedDot =
            "digraph G {"
            "rankdir=\"LR\";"
            "0[shape=Mrecord, label=\"{ ^ | 0 }\"];"
            "1[shape=Mrecord, label=\"{ $ | 0 }\"];"
            "2[shape=Mrecord, label=\"{ G | 1 }\"];"
            "3[shape=Mrecord, label=\"{ A | 2 }\"];"
            "4[shape=Mrecord, label=\"{ T | 2 }\"];"
            "0->2 ;"
            "2->3 ;"
            "3->4 ;"
            "4->1 ;"
            "0->3 ;"
            "}";
        PoaConsensusTests::plotConsensus(pc, "deletion-at-beginning");
        EXPECT_EQ(expectedDot, erase_all_copy(dot, "\n"));
        EXPECT_EQ("AT", pc->Sequence);
        delete pc;
    }

    // Deletion in middle
    {
        vector<string> reads;
        reads += "GAT", "GT";
        const PacBio::Poa::PoaConsensus* pc = PoaConsensus::FindConsensus(reads, AlignMode::GLOBAL);
        string dot = pc->Graph.ToGraphViz();
        string expectedDot =
            "digraph G {"
            "rankdir=\"LR\";"
            "0[shape=Mrecord, label=\"{ ^ | 0 }\"];"
            "1[shape=Mrecord, label=\"{ $ | 0 }\"];"
            "2[shape=Mrecord, label=\"{ G | 2 }\"];"
            "3[shape=Mrecord, label=\"{ A | 1 }\"];"
            "4[shape=Mrecord, label=\"{ T | 2 }\"];"
            "0->2 ;"
            "2->3 ;"
            "3->4 ;"
            "4->1 ;"
            "2->4 ;"
            "}";
        PoaConsensusTests::plotConsensus(pc, "deletion-in-middle");
        EXPECT_EQ(expectedDot, erase_all_copy(dot, "\n"));
        delete pc;
    }

    // Deletion at end
    {
        vector<string> reads;
        reads += "GAT", "GA";
        const PacBio::Poa::PoaConsensus* pc = PoaConsensus::FindConsensus(reads, AlignMode::GLOBAL);
        string dot = pc->Graph.ToGraphViz();
        string expectedDot =
            "digraph G {"
            "rankdir=\"LR\";"
            "0[shape=Mrecord, label=\"{ ^ | 0 }\"];"
            "1[shape=Mrecord, label=\"{ $ | 0 }\"];"
            "2[shape=Mrecord, label=\"{ G | 2 }\"];"
            "3[shape=Mrecord, label=\"{ A | 2 }\"];"
            "4[shape=Mrecord, label=\"{ T | 1 }\"];"
            "0->2 ;"
            "2->3 ;"
            "3->4 ;"
            "4->1 ;"
            "3->1 ;"
            "}";
        PoaConsensusTests::plotConsensus(pc, "deletion-at-end");
        EXPECT_EQ(expectedDot, erase_all_copy(dot, "\n"));
        EXPECT_EQ("GA", pc->Sequence);
        delete pc;
    }
}

TEST(PoaConsensus, TestSimple)
{
    vector<string> reads;
    reads += "TTTACAGGATAGTCCAGT", "ACAGGATACCCCGTCCAGT", "ACAGGATAGTCCAGT",
        "TTTACAGGATAGTCCAGTCCCC", "TTTACAGGATTAGTCCAGT", "TTTACAGGATTAGGTCCCAGT",
        "TTTACAGGATAGTCCAGT";
    const PacBio::Poa::PoaConsensus* pc = PoaConsensus::FindConsensus(reads, AlignMode::GLOBAL);
    PoaConsensusTests::plotConsensus(pc, "simple");
    EXPECT_EQ("TTTACAGGATAGTCCAGT", pc->Sequence);
    delete pc;
}

TEST(PoaConsensus, TestOverhangSecond)
{
    vector<string> reads;
    reads += "TTTACAGGATAGTCCAGT", "TTTACAGGATAGTCCAGTAAA", "TTTACAGGATAGTCCAGTAAA";
    const PacBio::Poa::PoaConsensus* pc = PoaConsensus::FindConsensus(reads, AlignMode::GLOBAL);
    EXPECT_EQ("TTTACAGGATAGTCCAGTAAA", pc->Sequence);
    delete pc;
}

TEST(PoaConsensus, SmallSemiglobalTest)
{
    vector<string> reads;
    reads += "GGTGG", "GGTGG", "T";
    const PacBio::Poa::PoaConsensus* pc = PoaConsensus::FindConsensus(reads, AlignMode::SEMIGLOBAL);
    PoaConsensusTests::plotConsensus(pc, "small-semiglobal");
    string expectedDot =
        "digraph G {"
        "rankdir=\"LR\";"
        "0[shape=Mrecord, label=\"{ ^ | 0 }\"];"
        "1[shape=Mrecord, label=\"{ $ | 0 }\"];"
        "2[shape=Mrecord, label=\"{ G | 2 }\"];"
        "3[shape=Mrecord, label=\"{ G | 2 }\"];"
        "4[shape=Mrecord, label=\"{ T | 3 }\"];"
        "5[shape=Mrecord, label=\"{ G | 2 }\"];"
        "6[shape=Mrecord, label=\"{ G | 2 }\"];"
        "0->2 ;"
        "2->3 ;"
        "3->4 ;"
        "4->5 ;"
        "5->6 ;"
        "6->1 ;"
        "4->1 ;"
        "0->4 ;"
        "}";
    string dot = pc->Graph.ToGraphViz();
    EXPECT_EQ(expectedDot, erase_all_copy(dot, "\n"));
    EXPECT_EQ("GGTGG", pc->Sequence);
    delete pc;
}

TEST(PoaConsensus, SmallTilingTest)
{
    vector<string> reads;
    reads += "GGGGAAAA", "AAAATTTT", "TTTTCCCC", "CCCCAGGA";
    const PacBio::Poa::PoaConsensus* pc = PoaConsensus::FindConsensus(reads, AlignMode::SEMIGLOBAL);
    PoaConsensusTests::plotConsensus(pc, "small-tiling");
    EXPECT_EQ("GGGGAAAATTTTCCCCAGGA", pc->Sequence);
    delete pc;
}

TEST(PoaConsensus, TestVerboseGraphVizOutput)
{
    vector<string> reads;
    reads += "GGG", "TGGG";
    const PacBio::Poa::PoaConsensus* pc = PoaConsensus::FindConsensus(reads, AlignMode::GLOBAL);
    string dot = pc->Graph.ToGraphViz(PoaGraph::COLOR_NODES | PoaGraph::VERBOSE_NODES, pc);

    string expectedDot =
        "digraph G {"
        "rankdir=\"LR\";"
        "0[shape=Mrecord, label=\"{ { 0 | ^ } | { 0 | 0 } | { 0.00 | 0.00 } "
        "}\"];"
        "1[shape=Mrecord, label=\"{ { 1 | $ } | { 0 | 0 } | { 0.00 | 0.00 } "
        "}\"];"
        "2[shape=Mrecord, style=\"filled\", fillcolor=\"lightblue\" ,"
        " label=\"{ { 2 | G } | { 2 | 2 } | { 2.00 | 2.00 } }\"];"
        "3[shape=Mrecord, style=\"filled\", fillcolor=\"lightblue\" ,"
        " label=\"{ { 3 | G } | { 2 | 2 } | { 2.00 | 4.00 } }\"];"
        "4[shape=Mrecord, style=\"filled\", fillcolor=\"lightblue\" ,"
        " label=\"{ { 4 | G } | { 2 | 2 } | { 2.00 | 6.00 } }\"];"
        "5[shape=Mrecord, label=\"{ { 5 | T } | { 1 | 2 } | { -0.00 | -0.00 } "
        "}\"];"
        "0->2 ;"
        "2->3 ;"
        "3->4 ;"
        "4->1 ;"
        "5->2 ;"
        "0->5 ;}";

    EXPECT_EQ(expectedDot, erase_all_copy(dot, "\n"));
    delete pc;
}

TEST(PoaConsensus, TestLocalStaggered)
{
    // Adapted from Pat's C# test
    vector<string> reads;
    reads += "TTTACAGGATAGTGCCGCCAATCTTCCAGT", "GATACCCCGTGCCGCCAATCTTCCAGTATATACAGCACGAGTAGC",
        "ATAGTGCCGCCAATCTTCCAGTATATACAGCACGGAGTAGCATCACGTACGTACGTCTACACGTAATT",
        "ACGTCTACACGTAATTTTGGAGAGCCCTCTCTCACG", "ACACGTAATTTTGGAGAGCCCTCTCTTCACG",
        "AGGATAGTGCCGCCAATCTTCCAGTAATATACAGCACGGAGTAGCATCACGTACG",
        "ATAGTGCCGCCAATCTTCCAGTATATACAGCACGGAGTAGCATCACGTACGTACGTCTACACGT";

    // 4 is a magic number here.  the .NET code had an entrenched
    // assumption that sequences in POA were subreads from a ZMW, so
    // the minCoverage was (numReads - 3), under assumption that basal
    // coverage for CCS is (numReads-2) (beginning, end read).
    // Application has to provide a sensible minCoverage.
    const PacBio::Poa::PoaConsensus* pc = PoaConsensus::FindConsensus(reads, AlignMode::LOCAL, 4);
    PoaConsensusTests::plotConsensus(pc, "local-staggered", false);
    EXPECT_EQ("ATAGTGCCGCCAATCTTCCAGTATATACAGCACGGAGTAGCATCACGTACGTACGTCTACACGTAATT", pc->Sequence);
    delete pc;
}

TEST(PoaConsensus, TestLongInsert)
{
    // Adapted from Pat's C# test
    vector<string> reads;
    reads +=
        "TTTACAGGATAGTGCCGCCAATCTTCCAGTGATACCCCGTGCCGCCAATCTTCCAGTATATACAGCACGA"
        "GGTAGC",
        "TTTACAGGATAGTGCCGGCCAATCTTCCAGTGATACCCCGTGCCGCCAATCTTCCAGTATATACAGCACG"
        "AGTAGC",
        "TTGTACAGGATAGTGCCGCCAATCTTCCAGTGATGGGGGGGGGGGGGGGGGGGGGGGGGGGACCCCGTGC"
        "CGCCAATCTTCCAGTATATACAGCACGAGTAGC";
    const PacBio::Poa::PoaConsensus* pc = PoaConsensus::FindConsensus(reads, AlignMode::GLOBAL);
    EXPECT_EQ(
        "TTTACAGGATAGTGCCGCCAATCTTCCAGTGATACCCCGTGCCGCCAATCTTCCAGTATATACAGCACGA"
        "GTAGC",
        pc->Sequence);
    delete pc;
}

TEST(PoaConsensus, TestSpanningReads)
{
    string read1 = "GAAAG";
    string read2 = "GATAG";
    vector<string> reads{read1, read1, read1, read2, read2, read2};
    const PacBio::Poa::PoaConsensus* pc = PoaConsensus::FindConsensus(reads, AlignMode::LOCAL);
    PoaConsensusTests::plotConsensus(pc, "spanning-reads");

    string dot = pc->Graph.ToGraphViz(PoaGraph::VERBOSE_NODES | PoaGraph::COLOR_NODES, pc);
    // We expect to get spanning reads of 6 for the middle A/T nodes,
    // but each only has 3 reads passing through.
    // The PoaGraph doesn't really expose an API, we can only check it
    // by looking at the GraphViz output.

    // clang-format off
    string expectedDot =
        "digraph G {"
        "rankdir=\"LR\";"
        "0[shape=Mrecord, label=\"{ { 0 | ^ } | { 0 | 0 } | { 0.00 | 0.00 } }\"];"
        "1[shape=Mrecord, label=\"{ { 1 | $ } | { 0 | 0 } | { 0.00 | 0.00 } }\"];"
        "2[shape=Mrecord, style=\"filled\", fillcolor=\"lightblue\" , label=\"{ { 2 | G } | { 6 | 6 } | { 6.00 | 6.00 } }\"];"
        "3[shape=Mrecord, style=\"filled\", fillcolor=\"lightblue\" , label=\"{ { 3 | A } | { 6 | 6 } | { 6.00 | 12.00 } }\"];"
        "4[shape=Mrecord, style=\"filled\", fillcolor=\"lightblue\" , label=\"{ { 4 | A } | { 3 | 6 } | { -0.00 | 12.00 } }\"];"
        "5[shape=Mrecord, style=\"filled\", fillcolor=\"lightblue\" , label=\"{ { 5 | A } | { 6 | 6 } | { 6.00 | 18.00 } }\"];"
        "6[shape=Mrecord, style=\"filled\", fillcolor=\"lightblue\" , label=\"{ { 6 | G } | { 6 | 6 } | { 6.00 | 24.00 } }\"];"
        "7[shape=Mrecord, label=\"{ { 7 | T } | { 3 | 6 } | { -0.00 | 12.00 } }\"];"
        "0->2 ;"
        "2->3 ;"
        "3->4 ;"
        "4->5 ;"
        "5->6 ;"
        "6->1 ;"
        "7->5 ;"
        "3->7 ;"
        "}";
    // clang-format on

    EXPECT_EQ(expectedDot, erase_all_copy(dot, "\n"));
    delete pc;
}

#if 0
TEST(PoaConsensus, TestMutations)
{
    using ::testing::ElementsAreArray;

    vector<string> reads;
    reads += "TGATTACAT",
             "TGATTACAT",
             "TGATTCAT",    // Deletion @ 5
             "TGATTATAT",   // Substitution @ 6
             "TGATTGACAT";  // Insertion @ 5

    const PacBio::Poa::PoaConsensus* pc = PoaConsensus::FindConsensus(reads, GLOBAL);

    const vector<ScoredMutation>* scoredMutations = pc->Mutations();
    vector<string> variantDescriptions;
    foreach (const ScoredMutation& scoredMutation, *scoredMutations)
    {
        variantDescriptions.push_back(scoredMutation.ToString());
    }
    sort(variantDescriptions.begin(), variantDescriptions.end());
    const char* expectedDescriptions[] = { "Deletion @5:6 -3.00",
                                           "Insertion (G) @5 -3.00",
                                           "Substitution (T) @6:7 -3.00" };
    ASSERT_THAT(variantDescriptions, ElementsAreArray(expectedDescriptions));
    delete pc;
}
#endif

TEST(PoaConsensus, NondeterminismRegressionTest)
{
    //
    // This is a regression test for a real-world case of
    // nondeterminism found in the POA on a quiver job on Staph.
    //
    vector<string> reads;
    reads +=
        "TATCAATCAACGAAATTCGCCAATTCCGTCATGAATGTCAATATCTAACTACACTTTAGAATACATTCTT"
        "TGACATGCCTGGCCTATTGATATTTCAATAAAATCAGACTATAAAGACAACTTACAAATGATCCTATAAA"
        "TTAAAGATCGAGAATCTAAAGAGTGAAATTAAAGCTAATTACTGCTTTAAAAATTTTACGTGCACACAAA"
        "AATGAATTTATCCTCATTATATCGAAAATACCATGAAGTATAGTAAGCTAACTTGAATATGATCATTAAT"
        "CGGCTATATGATTATTTTGATAATGCAATGAGCATCAATCTGAATTTATGACCTATCATTCGCGTTGCAT"
        "TTATTGAAGTGAAAATTCATGTACGCTTTTTTATTTTATTAATATAATCCTTGATATTGGTTATATACCA"
        "CGCTGTCACATAATTTTCAATAAATTTTTCTACTAAATGAAGTGTCTGTTATCTATCAC";
    reads +=
        "TATCAACAACGAAAATGCGCAGTTACGTCATGATTTATGTCAAATAATCTAAACGACACTTTCAGAAATA"
        "AATACATTCGAGAAGATGAATGCCTGGCGCAAAGTGATTATTTCAATAAAATATTTGTACCTTGAAAGAC"
        "AATTTACAAATGAATGCTATAAAATTTAAATGGATCCGGAGAATCTTTAAAGTACGTGAAATTAAAGGCT"
        "AAGATTACTGCGAAAAATTTTCGTGCACAAGAAATGAATGTTCCAGATTAGTATCGGAAAATAAGCCATG"
        "AAGAAGCTAGCATTAACTTGAATATGATCGATTTAATCGGCAGTATTGGTAATTATCTTGATAAGCAATT"
        "GAGCATCAACTGAAATTGAATGACTCTACATGCCTCGCTGAGTATGCGATTTATTGAAAGTGAAATTCAG"
        "TAAAGTTTATTGTTATGAATAAATGCGTACTTGGATGAATATCCCGACGGTAGTTCAAGTGTAAATGGAG"
        "TGAGGGGGTTCTTTCTTATAGAATAGTTTTATACTACTGATAAGGTGTAACCTGAGTGAGTCGTGATTTT"
        "AGAGTTACTTGCGAAC";

    std::set<string> answers;
    for (int run = 0; run < 100; run++) {
        const PacBio::Poa::PoaConsensus* pc = PoaConsensus::FindConsensus(reads, AlignMode::GLOBAL);
#if 0
        char fname[100];
        sprintf(fname, "/tmp/gr%03d.dot", run);
        pc->WriteGraphVizFile(fname, (PoaGraph::VERBOSE_NODES | PoaGraph::COLOR_NODES));
#endif
        answers.insert(pc->Sequence);
        delete pc;
    }
    ASSERT_EQ(1, answers.size());
}
