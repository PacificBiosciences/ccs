// Author: Derek Barnett

#include <gtest/gtest.h>

#include <algorithm>
#include <string>
#include <vector>

#include <pacbio/genomicconsensus/plurality/Plurality.h>

TEST(GenomicConsensusTests, aligned_ref_from_cigar)
{
    using Plurality = PacBio::GenomicConsensus::Plurality;

    {  // all match
        const std::string cigarStr = "5=";
        const auto cigar = PacBio::BAM::Cigar::FromStdString(cigarStr);
        const std::string inputRef = "ACGTA";
        const std::string expectedRef = "ACGTA";

        const auto alignedRef = Plurality::AlignedReference(cigar, inputRef);
        EXPECT_EQ(expectedRef, alignedRef);
    }
    {  // contains insertions
        const std::string cigarStr = "5=2I3=";
        const auto cigar = PacBio::BAM::Cigar::FromStdString(cigarStr);
        const std::string inputRef = "ACGTACCC";
        const std::string expectedRef = "ACGTA--CCC";

        const auto alignedRef = Plurality::AlignedReference(cigar, inputRef);
        EXPECT_EQ(expectedRef, alignedRef);
    }
}

TEST(GenomicConsensusTests, basecalls_from_alignment)
{
    using Plurality = PacBio::GenomicConsensus::Plurality;

    const std::string alnRef = "AC--GTAA-T";
    const std::string alnRead = "ACGGGT-TTT";
    const std::vector<std::string> expected{"A", "C", "GGG", "T", "-", "T", "TT"};

    const auto basecalls = Plurality::BaseCallsForAlignment(alnRead, alnRef, 10);
    EXPECT_TRUE(std::equal(expected.cbegin(), expected.cend(), basecalls.begin()));
}

TEST(GenomicConsensusTests, top_alleles_from_matrix_one_read)
{
    using Plurality = PacBio::GenomicConsensus::Plurality;
    using BaseCallsMatrix = Plurality::BaseCallsMatrix;

    const size_t windowLength = 7;
    const BaseCallsMatrix matrix{{"A", "C", "GGG", "T", "-", "T", "TT"}};

    const auto top2 = Plurality::TopAllelesFromMatrix(matrix, windowLength);
    ASSERT_EQ(windowLength, top2.size());

    const auto& top2_at_0 = top2.at(0);
    EXPECT_EQ("A", top2_at_0.firstAllele.bases);
    EXPECT_EQ(1, top2_at_0.firstAllele.frequency);
    EXPECT_EQ("", top2_at_0.secondAllele.bases);
    EXPECT_EQ(0, top2_at_0.secondAllele.frequency);
    EXPECT_EQ(1, top2_at_0.totalCoverage);

    const auto& top2_at_2 = top2.at(2);
    EXPECT_EQ("GGG", top2_at_2.firstAllele.bases);
    EXPECT_EQ(1, top2_at_2.firstAllele.frequency);
    EXPECT_EQ("", top2_at_2.secondAllele.bases);
    EXPECT_EQ(0, top2_at_2.secondAllele.frequency);
    EXPECT_EQ(1, top2_at_2.totalCoverage);

    const auto& top2_at_4 = top2.at(4);
    EXPECT_EQ("-", top2_at_4.firstAllele.bases);
    EXPECT_EQ(1, top2_at_4.firstAllele.frequency);
    EXPECT_EQ("", top2_at_4.secondAllele.bases);
    EXPECT_EQ(0, top2_at_4.secondAllele.frequency);
    EXPECT_EQ(1, top2_at_4.totalCoverage);

    const auto& top2_at_6 = top2.at(6);
    EXPECT_EQ("TT", top2_at_6.firstAllele.bases);
    EXPECT_EQ(1, top2_at_6.firstAllele.frequency);
    EXPECT_EQ("", top2_at_6.secondAllele.bases);
    EXPECT_EQ(0, top2_at_6.secondAllele.frequency);
    EXPECT_EQ(1, top2_at_6.totalCoverage);
}

TEST(GenomicConsensusTests, top_alleles_from_matrix_multiple_reads)
{
    using Plurality = PacBio::GenomicConsensus::Plurality;
    using BaseCallsMatrix = Plurality::BaseCallsMatrix;

    //    ref:   AC--GTAA-T
    //    read0: ACGGGT-TTT
    //    read1: ACGGGTAA-T
    //    read2: AC-GG-AA-T

    const size_t windowLength = 7;
    const BaseCallsMatrix matrix{{"A", "C", "GGG", "T", "-", "T", "TT"},
                                 {"A", "C", "GGG", "T", "A", "A", "T"},
                                 {"A", "C", "GG", "-", "A", "A", "T"}};

    const auto top2 = Plurality::TopAllelesFromMatrix(matrix, windowLength);
    ASSERT_EQ(windowLength, top2.size());

    const auto& top2_at_0 = top2.at(0);
    EXPECT_EQ("A", top2_at_0.firstAllele.bases);
    EXPECT_EQ(3, top2_at_0.firstAllele.frequency);
    EXPECT_EQ("", top2_at_0.secondAllele.bases);
    EXPECT_EQ(0, top2_at_0.secondAllele.frequency);
    EXPECT_EQ(3, top2_at_0.totalCoverage);

    const auto& top2_at_2 = top2.at(2);
    EXPECT_EQ("GGG", top2_at_2.firstAllele.bases);
    EXPECT_EQ(2, top2_at_2.firstAllele.frequency);
    EXPECT_EQ("GG", top2_at_2.secondAllele.bases);
    EXPECT_EQ(1, top2_at_2.secondAllele.frequency);
    EXPECT_EQ(3, top2_at_2.totalCoverage);

    const auto& top2_at_4 = top2.at(4);
    EXPECT_EQ("A", top2_at_4.firstAllele.bases);
    EXPECT_EQ(2, top2_at_4.firstAllele.frequency);
    EXPECT_EQ("-", top2_at_4.secondAllele.bases);
    EXPECT_EQ(1, top2_at_4.secondAllele.frequency);
    EXPECT_EQ(3, top2_at_4.totalCoverage);

    const auto& top2_at_6 = top2.at(6);
    EXPECT_EQ("T", top2_at_6.firstAllele.bases);
    EXPECT_EQ(2, top2_at_6.firstAllele.frequency);
    EXPECT_EQ("TT", top2_at_6.secondAllele.bases);
    EXPECT_EQ(1, top2_at_6.secondAllele.frequency);
    EXPECT_EQ(3, top2_at_6.totalCoverage);
}

TEST(GenomicConsensusTests, posterior_confidences)
{
    using Plurality = PacBio::GenomicConsensus::Plurality;

    const size_t depth = 20;
    const size_t cssFreq = 10;
    const size_t altFreq = 5;
    const double cap = 40.0;

    {  // haploid

        const bool diploid = false;
        const uint8_t expectedCssConf = 40;
        const uint8_t expectedHetConf = 0;

        const auto confidences =
            Plurality::PosteriorConfidences(depth, cssFreq, altFreq, diploid, cap);
        EXPECT_EQ(expectedCssConf, confidences.consensusConfidence);
        EXPECT_EQ(expectedHetConf, confidences.heterozygousConfidence);
    }
    {  // diploid

        // TODO: math is off on diploid code path

        //        const bool diploid = true;
        //        const uint8_t expectedCssConf = 40;
        //        const uint8_t expectedHetConf = 18;

        //        const auto confidences = Plurality::PosteriorConfidences(depth, cssFreq, altFreq, diploid, cap);
        //        EXPECT_EQ(expectedCssConf, confidences.consensusConfidence);
        //        EXPECT_EQ(expectedHetConf, confidences.heterozygousConfidence);
    }
}

TEST(GenomicConsensusTests, variants_from_ref_and_read)
{
    using Plurality = PacBio::GenomicConsensus::Plurality;
    using Allele = Plurality::Allele;

    // REF:   G
    // READ: AC
    //   => insertion(A), substitution(G->C)

    const std::string refName{"contig_foo"};
    const size_t refStart = 20;
    const char refBase = 'G';
    const Allele readAllele{"AC", 40};
    const uint8_t confidence = 40;
    const size_t coverage = 40;
    const char refPrev = 'T';
    const char readPrev = 'T';

    const auto variants = Plurality::VariantsFromRefAndRead(
        refName, refStart, refBase, readAllele, confidence, coverage, refPrev, readPrev);

    ASSERT_EQ(2, variants.size());

    const auto& v1 = variants.at(0);
    EXPECT_EQ(refName, v1.refName);
    EXPECT_EQ(refStart, v1.refStart);
    EXPECT_EQ("", v1.refSeq);
    EXPECT_EQ("A", v1.readSeq1);
    EXPECT_EQ(40, v1.frequency1.get());
    EXPECT_EQ(confidence, v1.confidence.get());
    EXPECT_EQ(coverage, v1.coverage.get());
    EXPECT_EQ(refPrev, v1.refPrev);
    EXPECT_EQ(readPrev, v1.readPrev);

    const auto& v2 = variants.at(1);
    EXPECT_EQ(refName, v2.refName);
    EXPECT_EQ(refStart, v2.refStart);
    EXPECT_EQ("G", v2.refSeq);
    EXPECT_EQ("C", v2.readSeq1);
    EXPECT_EQ(40, v2.frequency1.get());
    EXPECT_EQ(confidence, v2.confidence.get());
    EXPECT_EQ(coverage, v2.coverage.get());
    EXPECT_EQ(refPrev, v2.refPrev);
    EXPECT_EQ('T', v2.readPrev);
}

TEST(GenomicConsensusTests, variants_from_ref_and_reads)
{
    using Plurality = PacBio::GenomicConsensus::Plurality;
    using Allele = Plurality::Allele;

    // REF:  G
    // CSS: AC
    // ALT: AT
    //   => insertion(A)|insertion(A), substitution(G->C)|substitution(G->T)

    const std::string refName{"contig_foo"};
    const size_t refStart = 20;
    const char refBase = 'G';
    const Allele cssAllele{"AC", 40};
    const Allele altAllele{"AT", 30};
    const uint8_t confidence = 40;
    const size_t coverage = 70;
    const char refPrev = 'T';
    const char readPrev = 'T';

    const auto variants = Plurality::VariantsFromRefAndReads(
        refName, refStart, refBase, cssAllele, altAllele, confidence, coverage, refPrev, readPrev);

    ASSERT_EQ(2, variants.size());

    const auto& v1 = variants.at(0);
    EXPECT_EQ(refName, v1.refName);
    EXPECT_EQ(refStart, v1.refStart);
    EXPECT_EQ("", v1.refSeq);
    EXPECT_EQ("A", v1.readSeq1);
    EXPECT_EQ(40, v1.frequency1.get());
    EXPECT_EQ("A", v1.readSeq2.get());
    EXPECT_EQ(30, v1.frequency2.get());
    EXPECT_EQ(confidence, v1.confidence.get());
    EXPECT_EQ(coverage, v1.coverage.get());
    EXPECT_EQ(refPrev, v1.refPrev);
    EXPECT_EQ(readPrev, v1.readPrev);

    const auto& v2 = variants.at(1);
    EXPECT_EQ(refName, v2.refName);
    EXPECT_EQ(refStart, v2.refStart);
    EXPECT_EQ("G", v2.refSeq);
    EXPECT_EQ("C", v2.readSeq1);
    EXPECT_EQ(40, v2.frequency1.get());
    EXPECT_EQ("T", v2.readSeq2.get());
    EXPECT_EQ(30, v2.frequency2.get());
    EXPECT_EQ(confidence, v2.confidence.get());
    EXPECT_EQ(coverage, v2.coverage.get());
    EXPECT_EQ(refPrev, v2.refPrev);
    EXPECT_EQ('T', v2.readPrev);
}

TEST(GenomicConsensusTests, is_all_upper)
{
    using Plurality = PacBio::GenomicConsensus::Plurality;

    const std::string empty;
    const std::string single_lower = {"a"};
    const std::string single_upper{"A"};
    const std::string all_lower{"aaa"};
    const std::string all_upper{"AAA"};
    const std::string mixed{"AaA"};

    EXPECT_TRUE(Plurality::IsAllUpper(single_upper));
    EXPECT_TRUE(Plurality::IsAllUpper(all_upper));

    EXPECT_FALSE(Plurality::IsAllUpper(empty));
    EXPECT_FALSE(Plurality::IsAllUpper(single_lower));
    EXPECT_FALSE(Plurality::IsAllUpper(all_lower));
    EXPECT_FALSE(Plurality::IsAllUpper(mixed));
}
