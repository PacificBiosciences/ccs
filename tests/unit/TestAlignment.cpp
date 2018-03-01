// Authors: David Alexander, Lance Hepler

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <pacbio/align/AffineAlignment.h>
#include <pacbio/align/AlignConfig.h>
#include <pacbio/align/LinearAlignment.h>
#include <pacbio/align/LocalAlignment.h>
#include <pacbio/align/PairwiseAlignment.h>

// fwd declarations
namespace PacBio {
namespace Align {
namespace internal {
bool Rewrite2L(std::string* target, std::string* query, std::string* transcript, size_t i);
bool Rewrite3L(std::string* target, std::string* query, std::string* transcript, size_t i);
bool Rewrite2R(std::string* target, std::string* query, std::string* transcript, size_t i);
bool Rewrite3R(std::string* target, std::string* query, std::string* transcript, size_t i);
}
}
}

using namespace PacBio::Align;  // NOLINT
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

    PairwiseAlignment a2("GATTA-CA", "CA-TAACA");
    EXPECT_EQ("RMDMMIMM", a2.Transcript());
    EXPECT_FLOAT_EQ(5. / 8, a2.Accuracy());
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
    EXPECT_FLOAT_EQ(2. / 7, a->Accuracy());
    delete a;
}

TEST(PairwiseAlignmentTests, TargetPositionsInQueryTest)
{
    // MMM -> 0123
    {
        int expected[] = {0, 1, 2, 3};
        ASSERT_THAT(TargetToQueryPositions("MMM"), ElementsAreArray(expected));
    }

    // DMM -> 0012, MDM -> 0112, MMD -> 0122,
    {
        int expected1[] = {0, 0, 1, 2};
        int expected2[] = {0, 1, 1, 2};
        int expected3[] = {0, 1, 2, 2};
        ASSERT_THAT(TargetToQueryPositions("DMM"), ElementsAreArray(expected1));
        ASSERT_THAT(TargetToQueryPositions("MDM"), ElementsAreArray(expected2));
        ASSERT_THAT(TargetToQueryPositions("MMD"), ElementsAreArray(expected3));
    }

    // IMM -> 123, MIM -> 023, MMI -> 013,
    {
        int expected1[] = {1, 2, 3};
        int expected2[] = {0, 2, 3};
        int expected3[] = {0, 1, 3};
        ASSERT_THAT(TargetToQueryPositions("IMM"), ElementsAreArray(expected1));
        ASSERT_THAT(TargetToQueryPositions("MIM"), ElementsAreArray(expected2));
        ASSERT_THAT(TargetToQueryPositions("MMI"), ElementsAreArray(expected3));
    }

    // MRM, MDIM -> 0123
    // MIDM -> 0223
    {
        int expected1[] = {0, 1, 2, 3};
        int expected2[] = {0, 2, 2, 3};
        ASSERT_THAT(TargetToQueryPositions("MRM"), ElementsAreArray(expected1));
        ASSERT_THAT(TargetToQueryPositions("MDIM"), ElementsAreArray(expected1));
        ASSERT_THAT(TargetToQueryPositions("MIDM"), ElementsAreArray(expected2));
    }
}

// ---------------- Alignment justification tests ----------------------

TEST(PairwiseAlignmentTests, Rewriting)
{
    using namespace PacBio::Align::internal;

    std::string t, q, x;
    // Rewrite2L
    {
        t = "ACCT";
        q = "ACCT";
        x = "MMMM";
        EXPECT_FALSE(Rewrite2L(&t, &q, &x, 1));
        EXPECT_EQ("ACCT", t);
        EXPECT_EQ("ACCT", q);
        EXPECT_EQ("MMMM", x);

        t = "ACGT";
        q = "AC-T";
        x = "MMDM";
        EXPECT_FALSE(Rewrite2L(&t, &q, &x, 1));
        EXPECT_EQ("ACGT", t);
        EXPECT_EQ("AC-T", q);
        EXPECT_EQ("MMDM", x);

        t = "ACCT";
        q = "A-CT";
        x = "MDMM";
        EXPECT_FALSE(Rewrite2L(&t, &q, &x, 1));
        EXPECT_EQ("ACCT", t);
        EXPECT_EQ("A-CT", q);
        EXPECT_EQ("MDMM", x);

        t = "A-CT";
        q = "ACCT";
        x = "MIMM";
        EXPECT_FALSE(Rewrite2L(&t, &q, &x, 1));
        EXPECT_EQ("A-CT", t);
        EXPECT_EQ("ACCT", q);
        EXPECT_EQ("MIMM", x);

        t = "ACCT";
        q = "AC-T";
        x = "MMDM";
        EXPECT_TRUE(Rewrite2L(&t, &q, &x, 1));
        EXPECT_EQ("ACCT", t);
        EXPECT_EQ("A-CT", q);
        EXPECT_EQ("MDMM", x);

        t = "AC-T";
        q = "ACCT";
        x = "MMIM";
        EXPECT_TRUE(Rewrite2L(&t, &q, &x, 1));
        EXPECT_EQ("A-CT", t);
        EXPECT_EQ("ACCT", q);
        EXPECT_EQ("MIMM", x);
    }

    // Rewrite3L
    {
        t = "ACGCT";
        q = "ACGCT";
        x = "MMMMM";
        EXPECT_FALSE(Rewrite3L(&t, &q, &x, 1));
        EXPECT_EQ("ACGCT", t);
        EXPECT_EQ("ACGCT", q);
        EXPECT_EQ("MMMMM", x);

        t = "ACGGT";
        q = "AC--T";
        x = "MMDDM";
        EXPECT_FALSE(Rewrite3L(&t, &q, &x, 1));
        EXPECT_EQ("ACGGT", t);
        EXPECT_EQ("AC--T", q);
        EXPECT_EQ("MMDDM", x);

        t = "ACGCT";
        q = "A--CT";
        x = "MDDMM";
        EXPECT_FALSE(Rewrite3L(&t, &q, &x, 1));
        EXPECT_EQ("ACGCT", t);
        EXPECT_EQ("A--CT", q);
        EXPECT_EQ("MDDMM", x);

        t = "A--CT";
        q = "ACGCT";
        x = "MIIMM";
        EXPECT_FALSE(Rewrite3L(&t, &q, &x, 1));
        EXPECT_EQ("A--CT", t);
        EXPECT_EQ("ACGCT", q);
        EXPECT_EQ("MIIMM", x);

        t = "ACGCT";
        q = "AC--T";
        x = "MMDDM";
        EXPECT_TRUE(Rewrite3L(&t, &q, &x, 1));
        EXPECT_EQ("ACGCT", t);
        EXPECT_EQ("A--CT", q);
        EXPECT_EQ("MDDMM", x);

        t = "AC--T";
        q = "ACGCT";
        x = "MMIIM";
        EXPECT_TRUE(Rewrite3L(&t, &q, &x, 1));
        EXPECT_EQ("A--CT", t);
        EXPECT_EQ("ACGCT", q);
        EXPECT_EQ("MIIMM", x);
    }

    // Rewrite2R
    {
        t = "ACCT";
        q = "ACCT";
        x = "MMMM";
        EXPECT_FALSE(Rewrite2R(&t, &q, &x, 1));
        EXPECT_EQ("ACCT", t);
        EXPECT_EQ("ACCT", q);
        EXPECT_EQ("MMMM", x);

        t = "ACGT";
        q = "AC-T";
        x = "MMDM";
        EXPECT_FALSE(Rewrite2R(&t, &q, &x, 1));
        EXPECT_EQ("ACGT", t);
        EXPECT_EQ("AC-T", q);
        EXPECT_EQ("MMDM", x);

        t = "ACCT";
        q = "AC-T";
        x = "MMDM";
        EXPECT_FALSE(Rewrite2R(&t, &q, &x, 1));
        EXPECT_EQ("ACCT", t);
        EXPECT_EQ("AC-T", q);
        EXPECT_EQ("MMDM", x);

        t = "AC-T";
        q = "ACCT";
        x = "MMIM";
        EXPECT_FALSE(Rewrite2R(&t, &q, &x, 1));
        EXPECT_EQ("AC-T", t);
        EXPECT_EQ("ACCT", q);
        EXPECT_EQ("MMIM", x);

        t = "ACCT";
        q = "A-CT";
        x = "MDMM";
        EXPECT_TRUE(Rewrite2R(&t, &q, &x, 1));
        EXPECT_EQ("ACCT", t);
        EXPECT_EQ("AC-T", q);
        EXPECT_EQ("MMDM", x);

        t = "A-CT";
        q = "ACCT";
        x = "MIMM";
        EXPECT_TRUE(Rewrite2R(&t, &q, &x, 1));
        EXPECT_EQ("AC-T", t);
        EXPECT_EQ("ACCT", q);
        EXPECT_EQ("MMIM", x);
    }

    // Rewrite3R
    {
        t = "ACGCT";
        q = "ACGCT";
        x = "MMMMM";
        EXPECT_FALSE(Rewrite3R(&t, &q, &x, 1));
        EXPECT_EQ("ACGCT", t);
        EXPECT_EQ("ACGCT", q);
        EXPECT_EQ("MMMMM", x);

        t = "ACGGT";
        q = "AC--T";
        x = "MMDDM";
        EXPECT_FALSE(Rewrite3R(&t, &q, &x, 1));
        EXPECT_EQ("ACGGT", t);
        EXPECT_EQ("AC--T", q);
        EXPECT_EQ("MMDDM", x);

        t = "ACGCT";
        q = "AC--T";
        x = "MMDDM";
        EXPECT_FALSE(Rewrite3R(&t, &q, &x, 1));
        EXPECT_EQ("ACGCT", t);
        EXPECT_EQ("AC--T", q);
        EXPECT_EQ("MMDDM", x);

        t = "AC--T";
        q = "ACGCT";
        x = "MMIIM";
        EXPECT_FALSE(Rewrite3R(&t, &q, &x, 1));
        EXPECT_EQ("AC--T", t);
        EXPECT_EQ("ACGCT", q);
        EXPECT_EQ("MMIIM", x);

        t = "ACGCT";
        q = "A--CT";
        x = "MDDMM";
        EXPECT_TRUE(Rewrite3R(&t, &q, &x, 1));
        EXPECT_EQ("ACGCT", t);
        EXPECT_EQ("AC--T", q);
        EXPECT_EQ("MMDDM", x);

        t = "A--CT";
        q = "ACGCT";
        x = "MIIMM";
        EXPECT_TRUE(Rewrite3R(&t, &q, &x, 1));
        EXPECT_EQ("AC--T", t);
        EXPECT_EQ("ACGCT", q);
        EXPECT_EQ("MMIIM", x);
    }
}

TEST(PairwiseAlignmentTests, JustifyTest)
{
    // deletion
    {
        PairwiseAlignment a = PairwiseAlignment("AAAAAA", "AAA-AA");

        a.Justify(LRType::LEFT);
        EXPECT_EQ("AAAAAA", a.Target());
        EXPECT_EQ("-AAAAA", a.Query());
        EXPECT_EQ("DMMMMM", a.Transcript());

        a.Justify(LRType::RIGHT);
        EXPECT_EQ("AAAAAA", a.Target());
        EXPECT_EQ("AAAAA-", a.Query());
        EXPECT_EQ("MMMMMD", a.Transcript());
    }

    // insertion
    {
        PairwiseAlignment a = PairwiseAlignment("A-AAAA", "AAAAAA");

        a.Justify(LRType::LEFT);
        EXPECT_EQ("-AAAAA", a.Target());
        EXPECT_EQ("AAAAAA", a.Query());
        EXPECT_EQ("IMMMMM", a.Transcript());

        a.Justify(LRType::RIGHT);
        EXPECT_EQ("AAAAA-", a.Target());
        EXPECT_EQ("AAAAAA", a.Query());
        EXPECT_EQ("MMMMMI", a.Transcript());
    }

    // interruption in homopolymer
    {
        PairwiseAlignment a = PairwiseAlignment("GATTTACA", "GAGT-ACA");

        a.Justify(LRType::LEFT);
        EXPECT_EQ("GATTTACA", a.Target());
        EXPECT_EQ("GAG-TACA", a.Query());
        EXPECT_EQ("MMRDMMMM", a.Transcript());

        a.Justify(LRType::RIGHT);
        EXPECT_EQ("GATTTACA", a.Target());
        EXPECT_EQ("GAGT-ACA", a.Query());
        EXPECT_EQ("MMRMDMMM", a.Transcript());
    }

    // double bases, adjacent
    {
        PairwiseAlignment a = PairwiseAlignment("AAAAAA", "AAA--A");

        a.Justify(LRType::LEFT);
        EXPECT_EQ("AAAAAA", a.Target());
        EXPECT_EQ("--AAAA", a.Query());
        EXPECT_EQ("DDMMMM", a.Transcript());

        a.Justify(LRType::RIGHT);
        EXPECT_EQ("AAAAAA", a.Target());
        EXPECT_EQ("AAAA--", a.Query());
        EXPECT_EQ("MMMMDD", a.Transcript());
    }

    // double bases, separated
    {
        PairwiseAlignment a = PairwiseAlignment("AAAAAA", "A-AA-A");

        a.Justify(LRType::LEFT);
        EXPECT_EQ("AAAAAA", a.Target());
        EXPECT_EQ("--AAAA", a.Query());
        EXPECT_EQ("DDMMMM", a.Transcript());

        a.Justify(LRType::RIGHT);
        EXPECT_EQ("AAAAAA", a.Target());
        EXPECT_EQ("AAAA--", a.Query());
        EXPECT_EQ("MMMMDD", a.Transcript());
    }

    // intervening insertion
    {
        PairwiseAlignment a = PairwiseAlignment("A----A", "AATAAA");

        a.Justify(LRType::LEFT);
        EXPECT_EQ("----AA", a.Target());
        EXPECT_EQ("AATAAA", a.Query());
        EXPECT_EQ("IIIIMM", a.Transcript());

        a.Justify(LRType::RIGHT);
        EXPECT_EQ("AA----", a.Target());
        EXPECT_EQ("AATAAA", a.Query());
        EXPECT_EQ("MMIIII", a.Transcript());
    }

    // intervening match
    {
        PairwiseAlignment a = PairwiseAlignment("A-T--A", "AATAAA");

        a.Justify(LRType::LEFT);
        EXPECT_EQ("-AT--A", a.Target());
        EXPECT_EQ("AATAAA", a.Query());
        EXPECT_EQ("IMMIIM", a.Transcript());

        a.Justify(LRType::RIGHT);
        EXPECT_EQ("A-TA--", a.Target());
        EXPECT_EQ("AATAAA", a.Query());
        EXPECT_EQ("MIMMII", a.Transcript());
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
    const char* target =
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

    const char* query =
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

    const char* expectedAlignedTarget =
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
    AlignConfig config(params, AlignMode::GLOBAL);

    int score, peerScore;
    PairwiseAlignment *a, *peerAlignment;

    a = AlignLinear("GATTACA", "GATTACA", &score);
    EXPECT_EQ("GATTACA", a->Target());
    EXPECT_EQ("GATTACA", a->Query());
    EXPECT_EQ("MMMMMMM", a->Transcript());
    EXPECT_EQ(14, score);
    delete a;

    a = AlignLinear("CGAC", "GAAAACGAC", &score);
    EXPECT_EQ("-----CGAC", a->Target());
    EXPECT_EQ("GAAAACGAC", a->Query());
    EXPECT_EQ("IIIIIMMMM", a->Transcript());
    delete a;

    a = AlignLinear(
        "CATCAGGTAAGAAAGTACGATGCTACAGCTTGTGACTGGTGCGGCACTTTTGGCTGAGTTTCCTGTCCACCTCATGTATTCTGCCCTAAC"
        "GTCGGTCTTCACGCCATTACTAGACCGACAAAATGGAACCGGGGCCTTAAACCCCGTTCGAGGCGTAGCAAGGAGATAGGGTTATGAACT"
        "CTCCCAGTCAATATACCAACACATCGTGGGACGGATTGCAGAGCGAATCTATCCGCGCTCGCATAATTTAGTGTTGATC",
        "CATCAGGTAAGAAAAGTACGATGCTACAGCTGTGACTGGTGCGGCACTTTTTGGCCTGAGTTTCCTGTCCACTCATGTATTCTGGCCCTA"
        "ACTTAGGTCGGTCTTCACGCCATTTACTAGCACGAAAACGACAAAATTGGAAGCCGGGGCCTAAACACCCGTTCGAGGCGGTAGCAAGGA"
        "GATTAGGGTTATGAACTCTCCCAGTCAATGATACAAACAATCGTGGGACGGAATTGCAGAGCGAATCTATCCGCGCTCAAGCATAATTTA"
        "GTGTTGATC",
        &score);
    delete a;

    a = AlignLinear(
        "CATCAGGTAAGAAAGTACGATGCTACAGCTTGTGACTGGTGCGGCACTTTTGGCTGAGTTTCCTGTCCACCTCATGTATTCTGCCCTAAC"
        "GTCGGTCTTCACGCCATTACTAGACCGACAAAATGGAAGCCGGGGCCTTAAACCCCGTTCGAGGCGTAGCAAGGAGATAGGGTTATGAAC"
        "TCTCCCAGTCAATATACCAACACATCGTGGGACGGATTGCAGAGCGAATCTATCCGCGCTCGCATAATTTAGTGTTGATC",
        "CCCCGGGAATCTCTAGAATGCATCAGGTAAGAAAGTAACGATGCTTACACTTGTGACTGGTTGCGGCACTTTTGGTGAGTTTCCTGTCCA"
        "CTCAATGTATTCTGCCTAACGTCGTGTCTTCACGCCATTTACTAGACCGAGAAGGAAATTGGAAGGCCCGGGGGCCTTAAACGCCCGTTC"
        "GAGCGTAGCTAAGGAGATAGGGTTATGAACTCTCCCAGTCATATAGCCAACACATCGTGGAACGGAATTGCAGAGCGAATCTATCCGCTG"
        "CTCGCATAAATTTAGTGTTGATCCATAAAGCTTGCTGAGGACTAGTAGCTT",
        &score);
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
    EXPECT_FLOAT_EQ(2. / 7, a->Accuracy());
    delete a;

    const char* ref =
        "GTATTTTAAATAAAAACATTAAGTTATGACGAAGAAGAACGGAAACGCCTTAAACCGGAAAATTTTCATAAATAGCGAAAACCCGCGAGG"
        "TCGCCGCCC";
    const char* read =
        "GTATTTTAAATAAAAAAACATTATAGTTTAATGAACGAGAATGAACGGTAATACGCCTTTAAAGCCTGAAATATTTTTCCATAAATGTAA"
        "TTTCTGTATATAATCTCCGCGAGTGTCTGCCGCCC";

    a = AlignLinear(ref, read, &score);
    peerAlignment = Align(ref, read, &peerScore, config);
    EXPECT_EQ(score, peerScore);

    delete a;
    delete peerAlignment;
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

// ------------------ Local alignment tests ---------------------

TEST(LocalAlignmentTests, Simple)
{
    const std::string target = "CAGCCTTTCTGACCCGGAAATCAAAATAGGCACAACAAA";
    const std::string query = "CTGAGCCGGTAAATC";

    const auto a = LocalAlign(target, query);

    EXPECT_EQ(8, a.TargetBegin());
    EXPECT_EQ(21, a.TargetEnd());
    EXPECT_EQ(0, a.QueryBegin());
    EXPECT_EQ(14, a.QueryEnd());
    EXPECT_EQ(2, a.NumMismatches());
    EXPECT_EQ(21, a.Score());
}

// --------------- Semi-Global alignment tests ------------------

TEST(SemiGlobalAlignmentTests, Simple)
{
    const std::string target = "CAGCCTTTCTGACCCGGAAATCAAAATAGGCACAACAAA";
    const std::string query = "CTGAGCCGGTAAATC";
    const AlignConfig cfg(AlignParams::Default(), AlignMode::SEMIGLOBAL);

    const auto pa = Align(target, query, cfg);

    EXPECT_EQ(13, pa->Matches());
    EXPECT_EQ(2, pa->Errors());
    EXPECT_EQ(7, pa->ReferenceStart());
    EXPECT_EQ(21, pa->ReferenceEnd());
}
