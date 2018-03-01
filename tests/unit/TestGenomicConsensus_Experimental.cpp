// Author: Derek Barnett

#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <fstream>
#include <iostream>

#include <gtest/gtest.h>

#include <pbbam/EntireFileQuery.h>

#include <pacbio/data/Interval.h>
#include <pacbio/genomicconsensus/experimental/Consensus.h>
#include <pacbio/genomicconsensus/experimental/ConsensusMode.h>
#include <pacbio/genomicconsensus/experimental/ConsensusModelFactory.h>
#include <pacbio/genomicconsensus/experimental/Filters.h>
#include <pacbio/genomicconsensus/experimental/Input.h>
#include <pacbio/genomicconsensus/experimental/Intervals.h>
#include <pacbio/genomicconsensus/experimental/NoCallStyle.h>
#include <pacbio/genomicconsensus/experimental/ReferenceWindow.h>
#include <pacbio/genomicconsensus/experimental/Settings.h>
#include <pacbio/genomicconsensus/experimental/Sorting.h>
#include <pacbio/genomicconsensus/experimental/SortingStrategy.h>
#include <pacbio/genomicconsensus/experimental/WindowResult.h>
#include <pacbio/genomicconsensus/experimental/WorkChunk.h>
#include <pacbio/genomicconsensus/experimental/Workflow.h>
#include <pacbio/genomicconsensus/experimental/arrow/ArrowModel.h>
#include <pacbio/genomicconsensus/experimental/plurality/PluralityModel.h>
#include <pacbio/genomicconsensus/experimental/poa/PoaModel.h>

#include "TestData.h"

// clang-format off

using Interval = PacBio::Data::Interval;
using namespace PacBio::GenomicConsensus::experimental;

namespace GenomicConsensusExperimentalTests {

static const std::string All4merBam{
    tests::DataDir + "/genomicconsensus/all4mer/out.aligned_subreads.bam"};

static const std::string All4merFasta{
    tests::DataDir + "/genomicconsensus/all4mer/All4mer.V2.01_Insert.fa"};

static const std::string ChimeraFasta{
    tests::DataDir + "/chimera_minimal.fasta"};

std::vector<Variant> MakeFilteringTestVariants()
{
    Variant v1;
    v1.coverage = 0;
    v1.confidence = 0;

    Variant v2;
    v2.coverage = 5;
    v2.confidence = 0;

    Variant v3;
    v3.coverage = 0;
    v3.confidence = 40;

    Variant v4;
    v4.coverage = 4;
    v4.confidence = 70;

    Variant v5;
    v5.coverage = 6;
    v5.confidence = 30;

    return { v1, v2, v3, v4, v5 };
}

const std::vector<Variant>& FilteringTestVariants()
{
    static const std::vector<Variant> variants = MakeFilteringTestVariants();
    return variants;
}

std::vector<PacBio::BAM::BamRecord> MakeFilterSortTestReads()
{
    std::vector<PacBio::BAM::BamRecord> reads;
    PacBio::BAM::EntireFileQuery query{All4merBam};
    for (const auto& read : query)
        reads.push_back(read);
    return reads;
}

const std::vector<PacBio::BAM::BamRecord>& FilterSortTestReads()
{
    static const std::vector<PacBio::BAM::BamRecord> reads = MakeFilterSortTestReads();
    return reads;
}

}  // namespace GenomicConsensusExperimentalTests

// -----------------------
// Consensus
// -----------------------

TEST(GenomicConsensusExperimentalTest, no_call_consensus_with_no_call_style)
{
    const ReferenceWindow window{"foo", {0,10}};
    const std::string seq = "ACGTACGTAC";
    const NoCallStyle style = NoCallStyle::NO_CALL;

    const auto noCall = Consensus::NoCallConsensus(style, window, seq);

    const ReferenceWindow expectedWindow{"foo", {0,10}};
    const std::string expectedSeq = "NNNNNNNNNN";
    const std::vector<uint8_t> expectedConfidence(10, 0);

    EXPECT_EQ(expectedWindow, noCall.window);
    EXPECT_EQ(expectedSeq, noCall.sequence);
    EXPECT_EQ(expectedConfidence, noCall.confidence);
}

TEST(GenomicConsensusExperimentalTest, no_call_consensus_with_no_call_style_from_empty_input)
{
    const ReferenceWindow window{"foo", {0,0}};
    const std::string seq;
    const NoCallStyle style = NoCallStyle::NO_CALL;

    const auto noCall = Consensus::NoCallConsensus(style, window, seq);

    const ReferenceWindow expectedWindow{"foo", {0,0}};
    EXPECT_EQ(expectedWindow, noCall.window);
    EXPECT_TRUE(noCall.sequence.empty());
    EXPECT_TRUE(noCall.confidence.empty());
}

TEST(GenomicConsensusExperimentalTest, no_call_consensus_with_reference_call_style)
{
    const ReferenceWindow window{"foo", {0,10}};
    const std::string seq = "ACGTACGTAC";
    const NoCallStyle style = NoCallStyle::REFERENCE;

    const auto noCall = Consensus::NoCallConsensus(style, window, seq);

    const ReferenceWindow expectedWindow{"foo", {0,10}};
    const std::string expectedSeq = "ACGTACGTAC";
    const std::vector<uint8_t> expectedConfidence(10, 0);

    EXPECT_EQ(expectedWindow, noCall.window);
    EXPECT_EQ(expectedSeq, noCall.sequence);
    EXPECT_EQ(expectedConfidence, noCall.confidence);
}

TEST(GenomicConsensusExperimentalTest, no_call_consensus_with_reference_style_from_empty_input)
{
    const ReferenceWindow window{"foo", {0,0}};
    const std::string seq;
    const NoCallStyle style = NoCallStyle::REFERENCE;

    const auto noCall = Consensus::NoCallConsensus(style, window, seq);

    const ReferenceWindow expectedWindow{"foo", {0,0}};
    EXPECT_EQ(expectedWindow, noCall.window);
    EXPECT_TRUE(noCall.sequence.empty());
    EXPECT_TRUE(noCall.confidence.empty());
}

TEST(GenomicConsensusExperimentalTest, no_call_consensus_with_lowercase_reference_call_style)
{
    const ReferenceWindow window{"foo", {0,10}};
    const std::string seq = "ACGTACGTAC";
    const NoCallStyle style = NoCallStyle::LOWERCASE_REFERENCE;

    const auto noCall = Consensus::NoCallConsensus(style, window, seq);

    const ReferenceWindow expectedWindow{"foo", {0,10}};
    const std::string expectedSeq = "acgtacgtac";
    const std::vector<uint8_t> expectedConfidence(10, 0);

    EXPECT_EQ(expectedWindow, noCall.window);
    EXPECT_EQ(expectedSeq, noCall.sequence);
    EXPECT_EQ(expectedConfidence, noCall.confidence);
}

TEST(GenomicConsensusExperimentalTest, no_call_consensus_with_lowercase_reference_style_from_empty_input)
{
    const ReferenceWindow window{"foo", {0,0}};
    const std::string seq;
    const NoCallStyle style = NoCallStyle::LOWERCASE_REFERENCE;

    const auto noCall = Consensus::NoCallConsensus(style, window, seq);

    const ReferenceWindow expectedWindow{"foo", {0,0}};
    EXPECT_EQ(expectedWindow, noCall.window);
    EXPECT_TRUE(noCall.sequence.empty());
    EXPECT_TRUE(noCall.confidence.empty());
}

TEST(GenomicConsensusExperimentalTest, joining_empty_consensi_throws)
{
    const std::vector<Consensus> empty;
    EXPECT_THROW(
    {
        Consensus::Join(empty);
    },
    std::runtime_error);
}

TEST(GenomicConsensusExperimentalTest, can_join_consensus)
{
    const std::string seq{"ACGTACGTAC"};
    const std::vector<uint8_t> conf(10, 42);
    const Consensus left { ReferenceWindow{"foo", {0,10}}, seq, conf };
    const Consensus right { ReferenceWindow{"foo", {10,20}}, seq, conf };

    const auto joined = Consensus::Join({left, right});

    const ReferenceWindow joinedWindow{ "foo", {0,20}};
    const std::string joinedSeq{"ACGTACGTACACGTACGTAC"};
    const std::vector<uint8_t> joinedConf(20, 42);
    EXPECT_EQ(joinedWindow, joined.window);
    EXPECT_EQ(joinedSeq, joined.sequence);
    EXPECT_EQ(joinedConf, joined.confidence);
}

TEST(GenomicConsensusExperimentalTest, can_compare_consensus)
{
    const std::string seq{"ACGTACGTAC"};
    const std::vector<uint8_t> conf(10, 42);
    const Consensus left { ReferenceWindow{"foo", {0,10}}, seq, conf };
    const Consensus right { ReferenceWindow{"foo", {100,110}}, seq, conf  };

    EXPECT_LT(left, right);
}

TEST(GenomicConsensusExperimentalTest, factory_creates_expected_type_from_mode)
{
    {   // arrow
        const auto base = ConsensusModelFactory::Create(ConsensusMode::ARROW);
        const auto arrowModel = dynamic_cast<ArrowModel*>(base.get());
        EXPECT_NE(nullptr, arrowModel);
    }
    {   // plurality
        const auto base = ConsensusModelFactory::Create(ConsensusMode::PLURALITY);
        const auto pluralityModel = dynamic_cast<PluralityModel*>(base.get());
        EXPECT_NE(nullptr, pluralityModel);
    }
    {   // poa
        const auto base = ConsensusModelFactory::Create(ConsensusMode::POA);
        const auto poaModel = dynamic_cast<PoaModel*>(base.get());
        EXPECT_NE(nullptr, poaModel);
    }

    EXPECT_THROW(
    {
        ConsensusModelFactory::Create(static_cast<ConsensusMode>(44));
    },
    std::runtime_error);
}

// -----------------------
// Filters
// -----------------------

TEST(GenomicConsensusExperimentalTest, filtering_alignments_with_zeroed_criteria_returns_all)
{
    // No real criteria, so all reads should pass.

    const float readStumpinessThreshold = 0.0f;
    const float minHqRegionSnr = 0.0f;
    const float minReadScore = 0.0f;
    auto reads = GenomicConsensusExperimentalTests::FilterSortTestReads();

    FilterAlignments(&reads, readStumpinessThreshold, minHqRegionSnr, minReadScore);

    EXPECT_EQ(507, reads.size());
}

TEST(GenomicConsensusExperimentalTest, filtering_alignments_with_read_stumpiness_threshold)
{
    // This stumpiness threshold doesn't make sense, but we need at least some
    // reads to fail to check that the filter is working.

    const float readStumpinessThreshold = 1.1f;
    const float minHqRegionSnr = 0.0f;
    const float minReadScore = 0.0f;
    auto reads = GenomicConsensusExperimentalTests::FilterSortTestReads();

    FilterAlignments(&reads, readStumpinessThreshold, minHqRegionSnr, minReadScore);

    EXPECT_EQ(42, reads.size());
    for (const auto& read : reads)
    {
        const auto readLength = read.AlignedEnd() - read.AlignedStart();
        const auto refLength = read.ReferenceEnd() - read.ReferenceStart();
        EXPECT_GE(readLength, (refLength * readStumpinessThreshold));
    }
}

TEST(GenomicConsensusExperimentalTest, filtering_alignments_with_min_snr)
{
    const float readStumpinessThreshold = 0.0f;
    const float minHqRegionSnr = 12.0f;
    const float minReadScore = 0.0f;
    auto reads = GenomicConsensusExperimentalTests::FilterSortTestReads();

    FilterAlignments(&reads, readStumpinessThreshold, minHqRegionSnr, minReadScore);

    EXPECT_EQ(121, reads.size());
    for (const auto& read : reads)
    {
        const auto snr = read.SignalToNoise();
        const auto lowestSnr = *std::min_element(snr.cbegin(), snr.cend());
        EXPECT_GE(lowestSnr, minHqRegionSnr);
    }
}

TEST(GenomicConsensusExperimentalTest, filtering_alignments_with_min_read_score)
{
    const float readStumpinessThreshold = 0.0f;
    const float minHqRegionSnr = 0.0f;
    const float minReadScore = 0.88f;
    auto reads = GenomicConsensusExperimentalTests::FilterSortTestReads();

    FilterAlignments(&reads, readStumpinessThreshold, minHqRegionSnr, minReadScore);

    EXPECT_EQ(153, reads.size());
    for (const auto& read : reads)
       EXPECT_GE(read.ReadAccuracy(), minReadScore);
}

TEST(GenomicConsensusExperimentalTest, filtering_variants_with_zeroed_criteria_returns_all)
{
    // No real criteria, so all variants should pass.

    const size_t minCoverage = 0;
    const size_t minConfidence = 0;
    auto variants = GenomicConsensusExperimentalTests::FilteringTestVariants();

    FilterVariants(&variants, minCoverage, minConfidence);

    EXPECT_EQ(5, variants.size());
}

TEST(GenomicConsensusExperimentalTest, filtering_variants_with_min_coverage)
{
    const size_t minCoverage = 5;
    const size_t minConfidence = 0;
    auto variants = GenomicConsensusExperimentalTests::FilteringTestVariants();

    FilterVariants(&variants, minCoverage, minConfidence);

    EXPECT_EQ(2, variants.size());
    for (const auto& v : variants)
        EXPECT_GE(v.coverage, minCoverage);
}

TEST(GenomicConsensusExperimentalTest, filtering_variants_with_min_confidece)
{
    const size_t minCoverage = 0;
    const size_t minConfidence = 40;
    auto variants = GenomicConsensusExperimentalTests::FilteringTestVariants();

    FilterVariants(&variants, minCoverage, minConfidence);

    EXPECT_EQ(2, variants.size());
    for (const auto& v : variants)
        EXPECT_GE(v.confidence, minConfidence);
}

// -----------------------
// Input
// -----------------------

TEST(GenomicConsensusExperimentalTest, reads_from_full_ref_window_using_default_settings)
{
    Settings settings;
    settings.inputFilename = GenomicConsensusExperimentalTests::All4merBam;
    settings.referenceFilename = GenomicConsensusExperimentalTests::All4merFasta;

    const Input input{settings};
    const ReferenceWindow window
    {
        "All4mer.V2.01_Insert",
        {0,500}
    };

    const auto reads = input.ReadsInWindow(window);

    EXPECT_EQ(100, reads.size()); // Settings::Defaults::MaxCoverage = 100
}

TEST(GenomicConsensusExperimentalTest, reads_from_window_using_default_settings)
{
    Settings settings;
    settings.inputFilename = GenomicConsensusExperimentalTests::All4merBam;
    settings.referenceFilename = GenomicConsensusExperimentalTests::All4merFasta;

    const Input input{settings};
    const ReferenceWindow window
    {
        "All4mer.V2.01_Insert",
        {5,20}
    };

    const auto reads = input.ReadsInWindow(window);

    EXPECT_EQ(100, reads.size()); // Settings::Defaults::MaxCoverage = 100
}

TEST(GenomicConsensusExperimentalTest, reads_from_window_using_all_relaxed_settings)
{
    Settings settings;
    settings.inputFilename = GenomicConsensusExperimentalTests::All4merBam;
    settings.referenceFilename = GenomicConsensusExperimentalTests::All4merFasta;

    settings.maxCoverage = 600;
    settings.readStumpinessThreshold = 0.0f;
    settings.minHqRegionSnr = 0.0f;
    settings.minReadScore = 0.0f;
    settings.minMapQV = 0;

    const Input input{settings};
    const ReferenceWindow window { "All4mer.V2.01_Insert", {5,20} };


    const auto reads = input.ReadsInWindow(window);

    EXPECT_EQ(494, reads.size()); // 13 reads start after window (507 - 13)
}

TEST(GenomicConsensusExperimentalTest, reads_from_window_using_relaxed_max_coverage)
{
    Settings settings;
    settings.inputFilename = GenomicConsensusExperimentalTests::All4merBam;
    settings.referenceFilename = GenomicConsensusExperimentalTests::All4merFasta;

    settings.maxCoverage = 600;

    const Input input{settings};
    const ReferenceWindow window { "All4mer.V2.01_Insert", {5,20} };

    const auto reads = input.ReadsInWindow(window);

    EXPECT_EQ(494, reads.size());
}

TEST(GenomicConsensusExperimentalTest, reads_from_full_ref_window_using_relaxed_max_coverage)
{
    Settings settings;
    settings.inputFilename = GenomicConsensusExperimentalTests::All4merBam;
    settings.referenceFilename = GenomicConsensusExperimentalTests::All4merFasta;

    settings.maxCoverage = 600;

    const Input input{settings};
    const ReferenceWindow window { "All4mer.V2.01_Insert", {0,500} };

    const auto reads = input.ReadsInWindow(window);

    EXPECT_EQ(507, reads.size()); // all reads
}

TEST(GenomicConsensusExperimentalTest, reads_from_full_ref_window_using_strict_map_qv)
{
    Settings settings;
    settings.inputFilename = GenomicConsensusExperimentalTests::All4merBam;
    settings.referenceFilename = GenomicConsensusExperimentalTests::All4merFasta;

    settings.maxCoverage = 600;
    settings.minMapQV = 255;

    const Input input{settings};
    const ReferenceWindow window { "All4mer.V2.01_Insert", {0,500} };

    const auto reads = input.ReadsInWindow(window);

    EXPECT_EQ(0, reads.size()); // all reads have (MAPQ == 254)
}

TEST(GenomicConsensusExperimentalTest, reads_from_full_ref_window_using_strict_stumpiness)
{
    // This stumpiness threshold doesn't make sense, but we need at least some
    // reads to fail to check that the filter is working.

    Settings settings;
    settings.inputFilename = GenomicConsensusExperimentalTests::All4merBam;
    settings.referenceFilename = GenomicConsensusExperimentalTests::All4merFasta;

    settings.readStumpinessThreshold = 1.1f;

    const Input input{settings};
    const ReferenceWindow window { "All4mer.V2.01_Insert", {0,500} };

    const auto reads = input.ReadsInWindow(window);

    EXPECT_EQ(42, reads.size());
}

TEST(GenomicConsensusExperimentalTest, reads_from_full_ref_window_using_strict_snr)
{
    Settings settings;
    settings.inputFilename = GenomicConsensusExperimentalTests::All4merBam;
    settings.referenceFilename = GenomicConsensusExperimentalTests::All4merFasta;
    settings.maxCoverage = 600;

    settings.minHqRegionSnr = 7.0f;

    const Input input{settings};
    const ReferenceWindow window { "All4mer.V2.01_Insert", {0,500} };

    const auto reads = input.ReadsInWindow(window);

    EXPECT_EQ(382, reads.size());
}

TEST(GenomicConsensusExperimentalTest, reads_from_full_ref_window_using_strict_read_score)
{
    Settings settings;
    settings.inputFilename = GenomicConsensusExperimentalTests::All4merBam;
    settings.referenceFilename = GenomicConsensusExperimentalTests::All4merFasta;
    settings.maxCoverage = 600;

    settings.minReadScore = 1.0f;

    const Input input{settings};
    const ReferenceWindow window { "All4mer.V2.01_Insert", {0,500} };

    const auto reads = input.ReadsInWindow(window);

    EXPECT_EQ(0, reads.size()); // all reads have (MAPQ == 254)
}

TEST(GenomicConsensusExperimentalTest, reads_from_full_ref_window_using_all_criteria)
{
    Settings settings;
    settings.inputFilename = GenomicConsensusExperimentalTests::All4merBam;
    settings.referenceFilename = GenomicConsensusExperimentalTests::All4merFasta;

    settings.maxCoverage = 600;
    settings.readStumpinessThreshold = 0.5f;
    settings.minHqRegionSnr = 7.0f;
    settings.minReadScore = 0.8f;

    const Input input{settings};
    const ReferenceWindow window { "All4mer.V2.01_Insert", {5,20} };

    const auto reads = input.ReadsInWindow(window);

    EXPECT_EQ(371, reads.size()); // 13 reads start after window (507 - 13)
}

TEST(GenomicConsensusExperimentalTest, reference_sequence_from_window)
{
    Settings settings;
    settings.referenceFilename = GenomicConsensusExperimentalTests::ChimeraFasta;
    const Input input{settings};
    const ReferenceWindow window {
        "Barcode0--0_Cluster1_Phase1_NumReads297",
        {10,20}
    };

    const auto seq = input.ReferenceInWindow(window);

    EXPECT_EQ("TTGCAGAAAC", seq);
}

TEST(GenomicConsensusExperimentalTest, sequence_length_from_fasta)
{
    Settings settings;
    settings.referenceFilename = GenomicConsensusExperimentalTests::ChimeraFasta;
    const Input input{settings};
    const std::string name = "Barcode0--0_Cluster1_Phase1_NumReads297";

    const auto seqLength = input.SequenceLength(name);

    EXPECT_EQ(3152, seqLength);
}

// -----------------------
// Intervals
// -----------------------

TEST(GenomicConsensusExperimentalTest, coverage_intervals_from_intervals)
{
    const auto window = Interval{0,100};
    const auto intervals  = std::vector<Interval>
    {
        Interval{0, 10},
        Interval{5, 20},
        Interval{30, 50},
        Interval{40, 50},
        Interval{70, 80},
        Interval{75, 85},
        Interval{75, 90}
    };

    const auto covIntervals = CoverageIntervals(window, intervals);

    // [0, 20)   : 2
    // [20, 30)  : 0
    // [30, 50)  : 2
    // [50, 70)  : 0
    // [70, 90)  : 3
    // [90, 100) : 0

    ASSERT_EQ(6, covIntervals.size());

    EXPECT_EQ(Interval(0,20), covIntervals.at(0).interval);
    EXPECT_EQ(Interval(20,30), covIntervals.at(1).interval);
    EXPECT_EQ(Interval(30,50), covIntervals.at(2).interval);
    EXPECT_EQ(Interval(50,70), covIntervals.at(3).interval);
    EXPECT_EQ(Interval(70,90), covIntervals.at(4).interval);
    EXPECT_EQ(Interval(90,100), covIntervals.at(5).interval);

    EXPECT_EQ(2, covIntervals.at(0).coverage);
    EXPECT_EQ(0, covIntervals.at(1).coverage);
    EXPECT_EQ(2, covIntervals.at(2).coverage);
    EXPECT_EQ(0, covIntervals.at(3).coverage);
    EXPECT_EQ(3, covIntervals.at(4).coverage);
    EXPECT_EQ(0, covIntervals.at(5).coverage);
}

TEST(GenomicConsensusExperimentalTest, coverage_intervals_from_empty_input_intervals_is_window_with_zero_coverage)
{
    const auto window = Interval{0,100};
    const auto intervals  = std::vector<Interval>{};

    const auto covIntervals = CoverageIntervals(window, intervals);

    ASSERT_EQ(1, covIntervals.size());
    EXPECT_EQ(window, covIntervals.at(0).interval);
    EXPECT_EQ(0, covIntervals.at(0).coverage);
}

TEST(GenomicConsensusExperimentalTest, coverage_intervals_for_empty_window_throws)
{
    const auto window = Interval{};
    const auto intervals  = std::vector<Interval>
    {
        Interval{0, 10},
        Interval{5, 20},
        Interval{30, 50},
        Interval{40, 50},
        Interval{70, 80},
        Interval{75, 90}
    };

    EXPECT_THROW(
    {
        CoverageIntervals(window, intervals);
    },
    std::invalid_argument);
}

TEST(GenomicConsensusExperimentalTest, coverage_intervals_for_disjoint_window_throws)
{
    // window outside input range
    const auto window = Interval{200, 300};
    const auto intervals  = std::vector<Interval>{
        Interval{0, 10},
        Interval{5, 20},
        Interval{30, 50},
        Interval{40, 50},
        Interval{70, 80},
        Interval{75, 90}
    };

    EXPECT_THROW(
    {
        CoverageIntervals(window, intervals);
    },
    std::invalid_argument);
}

// ##
// FancyIntervals
// ##

TEST(GenomicConsensusExperimentalTest, fancy_intervals)
{
    const Interval window{0, 1000};
    const std::vector<Interval> readIntervals
    {
        Interval{0,400},
        Interval{100,600},
        Interval{200,800},
        Interval{200,500},
        Interval{300,700},
        Interval{450,550},
        Interval{600,1000},
        Interval{850,1000},
        Interval{850,1000},
        Interval{900, 1000},
        Interval{950,1000}
    };
    const size_t minCoverage = 5;

    const auto intervals = FancyIntervals(window, readIntervals, minCoverage);
    ASSERT_EQ(6, intervals.size());

    EXPECT_EQ(Interval(0,   300), intervals.at(0));     // hole
    EXPECT_EQ(Interval(300, 400), intervals.at(1));     // k-spanned
    EXPECT_EQ(Interval(400, 450), intervals.at(2));     // hole
    EXPECT_EQ(Interval(450, 500), intervals.at(3));     // k-spanned
    EXPECT_EQ(Interval(500, 950), intervals.at(4));     // hole
    EXPECT_EQ(Interval(950, 1000), intervals.at(5));    // k-spanned
}

TEST(GenomicConsensusExperimentalTest, all_read_intervals_from_empty_filter)
{
    const PacBio::BAM::BamFile bamFile{ GenomicConsensusExperimentalTests::All4merBam };
    const auto pbiFn = bamFile.PacBioIndexFilename();
    const PacBio::BAM::PbiRawData index{ pbiFn };
    const PacBio::BAM::PbiFilter filter;

    const auto intervals = FilteredIntervals(index, filter);

    EXPECT_EQ(507, intervals.size());
}

TEST(GenomicConsensusExperimentalTest, read_intervals_from_zmw_filter)
{
    const PacBio::BAM::BamFile bamFile{ GenomicConsensusExperimentalTests::All4merBam };
    const auto pbiFn = bamFile.PacBioIndexFilename();
    const PacBio::BAM::PbiRawData index{ pbiFn };
    const PacBio::BAM::PbiFilter filter
    {
        PacBio::BAM::PbiZmwFilter{28}
    };

    const auto intervals = FilteredIntervals(index, filter);

    EXPECT_EQ(11, intervals.size());
}

TEST(GenomicConsensusExperimentalTest, hole_in_empty_intervals_is_full_window)
{
    const Interval win{0,100};
    const std::vector<Interval> intervals;

    const auto holes = Holes(win, intervals);

    ASSERT_EQ(1, holes.size());
    EXPECT_EQ(win, holes.front());
}

TEST(GenomicConsensusExperimentalTest, no_holes_in_contiguous_intervals)
{
    const Interval win{0,100};
    const std::vector<Interval> intervals
    {
        Interval{0,50},
        Interval{50, 100}
    };

    const auto holes = Holes(win, intervals);

    EXPECT_TRUE(holes.empty());
}

TEST(GenomicConsensusExperimentalTest, holes_in_disjoint_intervals)
{
    const Interval win{0,100};
    const std::vector<Interval> intervals
    {
        Interval{10, 30},
        Interval{40, 60},
        Interval{70, 90}
    };

    const auto holes = Holes(win, intervals);

    ASSERT_EQ(4, holes.size());
    EXPECT_EQ(Interval(0, 10), holes.at(0));
    EXPECT_EQ(Interval(30, 40), holes.at(1));
    EXPECT_EQ(Interval(60, 70), holes.at(2));
    EXPECT_EQ(Interval(90, 100), holes.at(3));
}

TEST(GenomicConsensusExperimentalTest, kspanned_intervals_from_empty_read_intervals_is_empty_list)
{
    const Interval windowInterval {0, 1000};
    const std::vector<Interval> readIntervals{};
    const size_t minCoverage = 5;

    const auto intervals = KSpannedIntervals(windowInterval, readIntervals, minCoverage);

    EXPECT_EQ(0, intervals.size());
}

TEST(GenomicConsensusExperimentalTest, kspanned_intervals_from_empty_window_throws)
{
    const Interval windowInterval {};
    const std::vector<Interval> readIntervals
    {
        Interval{0,400},
        Interval{100,600},
        Interval{200,800},
        Interval{200,500},
        Interval{300,700},
        Interval{450,550},
        Interval{600,1000},
        Interval{850,1000},
        Interval{850,1000},
        Interval{900, 1000},
        Interval{950,1000}
    };
    const size_t minCoverage = 5;

    EXPECT_THROW(
    {
        auto intervals = KSpannedIntervals(windowInterval, readIntervals, minCoverage);
    }, std::invalid_argument);
}

TEST(GenomicConsensusExperimentalTest, kspanned_intervals_over_window)
{
    const Interval windowInterval {0, 1000};
    const std::vector<Interval> readIntervals
    {
        Interval{0,400},
        Interval{100,600},
        Interval{200,800},
        Interval{200,500},
        Interval{300,700},
        Interval{450,550},
        Interval{600,1000},
        Interval{850,1000},
        Interval{850,1000},
        Interval{900, 1000},
        Interval{950,1000}
    };
    const size_t minCoverage = 5;

    const auto intervals = KSpannedIntervals(windowInterval, readIntervals, minCoverage);
    ASSERT_EQ(3, intervals.size());

    EXPECT_EQ(Interval(300, 400), intervals.at(0));
    EXPECT_EQ(Interval(450, 500), intervals.at(1));
    EXPECT_EQ(Interval(950, 1000), intervals.at(2));
}

TEST(GenomicConsensusExperimentalTest, projecting_from_empty_intervals_is_window_with_zero_coverage)
{
    const auto window = ReferenceWindow{"", {0,100}};
    const auto intervals  = std::vector<Interval>{};

    const auto projection = ProjectIntoRange(intervals, window.interval);

    EXPECT_EQ(window.Length(), projection.size());
    for (const auto p : projection)
        EXPECT_EQ(0, p);
}

TEST(GenomicConsensusExperimentalTest, projecting_intervals_from_empty_window_is_empty_list)
{
    const auto window = ReferenceWindow{"", {}};
    const auto intervals  = std::vector<Interval>
    {
        Interval{0, 10},
        Interval{5, 20},
        Interval{30, 50},
        Interval{40, 50},
        Interval{70, 80},
        Interval{75, 90}
    };

    const auto projection = ProjectIntoRange(intervals, window.interval);
    EXPECT_TRUE(projection.empty());
}

TEST(GenomicConsensusExperimentalTest, projecting_intervals_to_window)
{
    const auto window = ReferenceWindow{"", {0,20}};
    const auto intervals  = std::vector<Interval>
    {
        Interval{2, 8},
        Interval{5, 7},
        Interval{10, 15},
        Interval{12, 17},
    };

    const std::vector<size_t> expected =
    {
        0,0,1,1,1,2,2,1,0,0,1,1,2,2,2,1,1,0,0,0
    };

    const auto projection = ProjectIntoRange(intervals, window.interval);
    ASSERT_EQ(expected.size(), projection.size());
    for (size_t i = 0; i < expected.size(); ++i)
        EXPECT_EQ(expected.at(i), projection.at(i));
}

TEST(GenomicConsensusExperimentalTest, splitting_intervals_yields_contiguous_intervals_of_span_size)
{
    const auto source = Interval{0, 100};
    const auto span = 20;

    const auto intervals = SplitInterval(source, span);

    ASSERT_EQ(5, intervals.size());
    EXPECT_EQ(Interval(0,20), intervals.at(0));
    EXPECT_EQ(Interval(20,40), intervals.at(1));
    EXPECT_EQ(Interval(40,60), intervals.at(2));
    EXPECT_EQ(Interval(60,80), intervals.at(3));
    EXPECT_EQ(Interval(80,100), intervals.at(4));
}

TEST(GenomicConsensusExperimentalTest, splitting_intervals_clips_to_bounds)
{
    const auto source = Interval{10, 100};
    const auto span = 20;

    const auto intervals = SplitInterval(source, span);

    ASSERT_EQ(5, intervals.size());
    EXPECT_EQ(Interval(10, 30), intervals.at(0));
    EXPECT_EQ(Interval(30, 50), intervals.at(1));
    EXPECT_EQ(Interval(50, 70), intervals.at(2));
    EXPECT_EQ(Interval(70, 90), intervals.at(3));
    EXPECT_EQ(Interval(90, 100), intervals.at(4));
}

TEST(GenomicConsensusExperimentalTest, splitting_intervals_on_empty_interval_returns_none)
{
    const auto source = Interval{};
    const auto span = 20;

    const auto intervals = SplitInterval(source, span);

    EXPECT_TRUE(intervals.empty());
}

TEST(GenomicConsensusExperimentalTest, splitting_intervals_with_span_too_small_returns_input_interval)
{
    const auto source = Interval{0, 5};
    const auto span = 20;

    const auto intervals = SplitInterval(source, span);

    ASSERT_EQ(1, intervals.size());
    EXPECT_EQ(Interval(0, 5), intervals.at(0));
}

// -----------------------
// IPoaModel - free fxns ??
// -----------------------

// -----------------------
// ReferenceWindow
// -----------------------

TEST(GenomicConsensusExperimentalTest, reference_windows_compare_equal)
{
    const ReferenceWindow window1{"foo", {0, 100}};
    const ReferenceWindow window2{"foo", {0, 100}};
    EXPECT_EQ(window1, window2);
}

TEST(GenomicConsensusExperimentalTest, reference_windows_compare_not_equal)
{
    {   // different name
        const ReferenceWindow window1{"foo", {0, 100}};
        const ReferenceWindow window2{"bar", {0, 100}};
        EXPECT_NE(window1, window2);
    }
    {   // different interval
        const ReferenceWindow window1{"foo", {0, 100}};
        const ReferenceWindow window2{"foo", {0, 90}};
        EXPECT_NE(window1, window2);
    }
}

TEST(GenomicConsensusExperimentalTest, reference_windows_compare_less_than)
{
    {   // name less-than
        const ReferenceWindow window1{"foo", {0, 100}};
        const ReferenceWindow window2{"bar", {0, 100}};
        EXPECT_LT(window2, window1);
    }
    {   // interval less-than
        const ReferenceWindow window1{"foo", {0, 100}};
        const ReferenceWindow window2{"foo", {0, 90}};
        EXPECT_NE(window2, window1);
    }
}

TEST(GenomicConsensusExperimentalTest, adjacent_reference_windows_are_contiguous)
{
    const ReferenceWindow left{"foo", {0, 100}};
    const ReferenceWindow right{"foo", {100, 200}};
    EXPECT_TRUE(AreContiguous({left,right}));
}

TEST(GenomicConsensusExperimentalTest, reference_windows_on_different_refs_are_not_contiguous)
{
    const ReferenceWindow left{"foo", {0, 100}};
    const ReferenceWindow right{"bar", {100, 200}};
    EXPECT_FALSE(AreContiguous({left,right}));
}

TEST(GenomicConsensusExperimentalTest, disjoint_reference_windows_are_not_contiguous)
{
    const ReferenceWindow left{"foo", {0, 100}};
    const ReferenceWindow right{"foo", {200, 300}};
    EXPECT_FALSE(AreContiguous({left,right}));
}

TEST(GenomicConsensusExperimentalTest, overlapping_reference_windows_are_not_contiguous)
{
    const ReferenceWindow left{"foo", {0, 100}};
    const ReferenceWindow right{"foo", {50, 150}};
    EXPECT_FALSE(AreContiguous({left,right}));
}

TEST(GenomicConsensusExperimentalTest, identical_reference_windows_are_not_contiguous)
{
    const ReferenceWindow left{"foo", {0, 100}};
    const ReferenceWindow right{"foo", {0, 200}};
    EXPECT_FALSE(AreContiguous({left,right}));
}

TEST(GenomicConsensusExperimentalTest, can_print_reference_window)
{
    const ReferenceWindow window{"foo", {0, 100}};
    std::stringstream s;
    s << window;

    const std::string expected{"foo [0, 100)"};
    EXPECT_EQ(expected, s.str());
}

// -----------------------
// Settings / CLI
// -----------------------

// ##
// CreateInterface
// ##

// ##
// ctor from CLI
// ##

// -----------------------
// Sorting
// -----------------------

TEST(GenomicConsensusExperimentalTest, sorted_reads_by_longest_and_strand_balanced)
{
    const auto& reads = GenomicConsensusExperimentalTests::FilterSortTestReads();
    const std::vector<PacBio::BAM::BamRecord> readsToSort{ reads.cbegin(), reads.cbegin()+10 };
    const ReferenceWindow window{"All4mer.V2.01_Insert", {0, 500}};

    const auto sortedReads =
        SortedReadsInWindow(readsToSort, window, SortingStrategy::LONGEST_AND_STRAND_BALANCED);

    ASSERT_TRUE(sortedReads.size() == 10);
    // 260
    EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/2409_2745",
              sortedReads.at(0).FullName());
    EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/1669_1990",
              sortedReads.at(1).FullName());
    EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/193_534",
              sortedReads.at(2).FullName());
    EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/3923_4231",
              sortedReads.at(3).FullName());
    EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/9763_10082",
              sortedReads.at(4).FullName());
    // 259
    EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/943_1260",
              sortedReads.at(5).FullName());
    EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/9022_9354",
              sortedReads.at(6).FullName());
    EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/10491_10819",
              sortedReads.at(7).FullName());
    // 258
    EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/3189_3513",
              sortedReads.at(8).FullName());
    EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/4643_4956",
              sortedReads.at(9).FullName());

    // TODO: check dataset where strand-balance differs from just "longest"
}

TEST(GenomicConsensusExperimentalTest, sorted_reads_by_longest)
{
    const auto& reads = GenomicConsensusExperimentalTests::FilterSortTestReads();
    const std::vector<PacBio::BAM::BamRecord> readsToSort{ reads.cbegin(), reads.cbegin()+10 };
    const ReferenceWindow window{"All4mer.V2.01_Insert", {0, 500}};

    const auto sortedReads =
            SortedReadsInWindow(readsToSort, window, SortingStrategy::LONGEST);

    ASSERT_TRUE(sortedReads.size() == 10);
    // 260
    EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/2409_2745",
              sortedReads.at(0).FullName());
    EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/1669_1990",
              sortedReads.at(1).FullName());
    EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/193_534",
              sortedReads.at(2).FullName());
    EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/3923_4231",
              sortedReads.at(3).FullName());
    EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/9763_10082",
              sortedReads.at(4).FullName());
    // 259
    EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/943_1260",
              sortedReads.at(5).FullName());
    EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/9022_9354",
              sortedReads.at(6).FullName());
    EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/10491_10819",
              sortedReads.at(7).FullName());
    // 258
    EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/3189_3513",
              sortedReads.at(8).FullName());
    EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/4643_4956",
              sortedReads.at(9).FullName());
}

TEST(GenomicConsensusExperimentalTest, sorted_reads_by_spanning)
{
    {   // all should span, so same as file_order
        const auto& reads = GenomicConsensusExperimentalTests::FilterSortTestReads();
        const std::vector<PacBio::BAM::BamRecord> readsToSort{ reads.cbegin(), reads.cbegin()+10 };
        const ReferenceWindow smallWindow{"All4mer.V2.01_Insert", {0, 250}};

        const auto sortedReads
                = SortedReadsInWindow(readsToSort, smallWindow, SortingStrategy::SPANNING);

        ASSERT_EQ(10, sortedReads.size());
        EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/2409_2745",
                  sortedReads.at(0).FullName());
        EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/1669_1990",
                  sortedReads.at(1).FullName());
        EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/193_534",
                  sortedReads.at(2).FullName());
        EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/3189_3513",
                  sortedReads.at(3).FullName());
        EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/3923_4231",
                  sortedReads.at(4).FullName());
        EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/4643_4956",
                  sortedReads.at(5).FullName());
        EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/943_1260",
                  sortedReads.at(6).FullName());
        EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/9022_9354",
                  sortedReads.at(7).FullName());
        EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/9763_10082",
                  sortedReads.at(8).FullName());
        EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/10491_10819",
                  sortedReads.at(9).FullName());
    }
    {   // none fully span, so empty result
        const auto& reads = GenomicConsensusExperimentalTests::FilterSortTestReads();
        const std::vector<PacBio::BAM::BamRecord> readsToSort{ reads.cbegin(), reads.cbegin()+10 };
        const ReferenceWindow window{"All4mer.V2.01_Insert", {0, 500}};

        const auto sortedReads
                = SortedReadsInWindow(readsToSort, window, SortingStrategy::SPANNING);

        EXPECT_TRUE(sortedReads.empty());
    }
}

TEST(GenomicConsensusExperimentalTest, sorted_reads_by_file_order)
{
    const auto& reads = GenomicConsensusExperimentalTests::FilterSortTestReads();
    const std::vector<PacBio::BAM::BamRecord> readsToSort{ reads.cbegin(), reads.cbegin()+10 };
    const ReferenceWindow window{"All4mer.V2.01_Insert", {0, 500}};

    const auto sortedReads = SortedReadsInWindow(readsToSort, window, SortingStrategy::FILE_ORDER);

    ASSERT_TRUE(sortedReads.size() == 10);
    EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/2409_2745",
              sortedReads.at(0).FullName());
    EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/1669_1990",
              sortedReads.at(1).FullName());
    EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/193_534",
              sortedReads.at(2).FullName());
    EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/3189_3513",
              sortedReads.at(3).FullName());
    EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/3923_4231",
              sortedReads.at(4).FullName());
    EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/4643_4956",
              sortedReads.at(5).FullName());
    EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/943_1260",
              sortedReads.at(6).FullName());
    EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/9022_9354",
              sortedReads.at(7).FullName());
    EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/9763_10082",
              sortedReads.at(8).FullName());
    EXPECT_EQ("m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/10491_10819",
              sortedReads.at(9).FullName());
}

// -----------------------
// Utils - refactor around? available elsewhere?
// -----------------------

// -----------------------
// Variant
// -----------------------

TEST(GenomicConsensusExperimentalTest, can_annotate_variant)
{
    Variant v;
    v.Annotate("key", "value");

    EXPECT_EQ("key",   v.annotations.begin()->first);
    EXPECT_EQ("value", v.annotations.begin()->second);
}

TEST(GenomicConsensusExperimentalTest, empty_alt_allele_on_variant_is_homozygous)
{
    Variant v;
    v.readSeq1 = "C";

    EXPECT_TRUE(v.IsHomozygous());
    EXPECT_FALSE(v.IsHeterozygous());
}

TEST(GenomicConsensusExperimentalTest, same_alt_allele_on_variant_is_homozygous)
{
    Variant v;
    v.readSeq1 = "C";
    v.readSeq2 = "C";

    EXPECT_TRUE(v.IsHomozygous());
    EXPECT_FALSE(v.IsHeterozygous());
}

TEST(GenomicConsensusExperimentalTest, alt_allele_on_variant_is_heterozygous)
{
    Variant v;
    v.readSeq1 = "C";
    v.readSeq2 = "G";

    EXPECT_TRUE(v.IsHeterozygous());
    EXPECT_FALSE(v.IsHomozygous());
}

TEST(GenomicConsensusExperimentalTest, variant_compare_ordering)
{
    {   // first by refName
        Variant lhs;
        lhs.refName = "aa";
        lhs.refStart = 3;
        lhs.refEnd   = 4;
        lhs.readSeq1 = "yy";

        Variant rhs;
        rhs.refName = "bb";
        rhs.refStart = 3;
        rhs.refEnd   = 4;
        rhs.readSeq1 = "yy";

        EXPECT_TRUE(lhs < rhs);
    }
    {   // then by refStart
        Variant lhs;
        lhs.refName = "aa";
        lhs.refStart = 2;
        lhs.refEnd   = 4;
        lhs.readSeq1 = "zz";

        Variant rhs;
        rhs.refName = "aa";
        rhs.refStart = 3;
        rhs.refEnd   = 4;
        rhs.readSeq1 = "yy";

        EXPECT_TRUE(lhs < rhs);
    }
    {   // then by refEnd
        Variant lhs;
        lhs.refName = "aa";
        lhs.refStart = 3;
        lhs.refEnd   = 4;
        lhs.readSeq1 = "zz";

        Variant rhs;
        rhs.refName = "aa";
        rhs.refStart = 3;
        rhs.refEnd   = 5;
        rhs.readSeq1 = "yy";

        EXPECT_TRUE(lhs < rhs);
    }
    {   // last, by readSeq1
        Variant lhs;
        lhs.refName = "aa";
        lhs.refStart = 3;
        lhs.refEnd   = 5;
        lhs.readSeq1 = "kk";

        Variant rhs;
        rhs.refName = "aa";
        rhs.refStart = 3;
        rhs.refEnd   = 5;
        rhs.readSeq1 = "pp";

        EXPECT_TRUE(lhs < rhs);
    }
}

// -----------------------
// Workflow
// -----------------------

TEST(GenomicConsensusExperimentalTest, enumerate_chunks_from_filter_windows)
{
    const std::string name = "foo";
    const size_t stride = 20;
    const std::vector<ReferenceWindow> windows =
    {
        ReferenceWindow{"foo", {0,100}},    // 5 chunks
        ReferenceWindow{"bar", {0,200}},
        ReferenceWindow{"baz", {300, 450}},
        ReferenceWindow{"foo", {700,800}},  // 5 chunks
        ReferenceWindow{"foo", {200, 400}}  // 10 chunks
    };

    const auto chunks = Workflow::EnumerateChunks(name, stride, windows);
    EXPECT_EQ(20, chunks.size());
}

TEST(GenomicConsensusExperimentalTest, enumerate_chunks_returns_none_from_empty_filter_windows)
{
    const std::string name = "foo";
    const size_t stride = 20;
    const std::vector<ReferenceWindow> windows;

    const auto chunks = Workflow::EnumerateChunks(name, stride, windows);
    EXPECT_TRUE(chunks.empty());
}

TEST(GenomicConsensusExperimentalTest, enumerate_windows_from_filter_windows)
{
    const std::string name = "foo";
    const std::vector<ReferenceWindow> filterWindows =
    {
        ReferenceWindow{"foo", {0,100}},
        ReferenceWindow{"bar", {0,200}},
        ReferenceWindow{"baz", {300, 450}},
        ReferenceWindow{"foo", {700,800}},
        ReferenceWindow{"foo", {200, 400}}
    };

    const auto windows = Workflow::EnumerateWindows(name, filterWindows);
    EXPECT_EQ(3, windows.size());
    for (const auto& win : windows)
        EXPECT_EQ("foo", win.name);
}

TEST(GenomicConsensusExperimentalTest, enumerate_windows_returns_none_from_empty_filter_windows)
{
    const std::string name = "foo";
    const std::vector<ReferenceWindow> filterWindows;

    const auto windows = Workflow::EnumerateWindows(name, filterWindows);
    EXPECT_TRUE(windows.empty());
}

TEST(GenomicConsensusExperimentalTest, enumerate_spans_returns_ref_from_settings_with_empty_filter_windows)
{
    Settings settings;
    settings.referenceFilename = GenomicConsensusExperimentalTests::ChimeraFasta;
    const std::string name = "Barcode0--0_Cluster1_Phase1_NumReads297";

    const auto windows = Workflow::EnumerateWindows(name, settings);
    ASSERT_EQ(1, windows.size());
    EXPECT_EQ(name, windows.at(0).name);
    EXPECT_EQ(0, windows.at(0).Start());
    EXPECT_EQ(3152, windows.at(0).End());
}

TEST(GenomicConsensusExperimentalTest, simple_chunks_from_ref_name)
{
    const std::string name = "Barcode0--0_Cluster1_Phase1_NumReads297";

    Settings settings;
    settings.referenceFilename = GenomicConsensusExperimentalTests::ChimeraFasta;
    settings.windowSpan = 100;

    const auto chunks = Workflow::SimpleChunks(name, settings);

    EXPECT_EQ(32, chunks.size());   // 3152 bp / 100 span
}

TEST(GenomicConsensusExperimentalTest, reference_names_from_file)
{
    Settings settings;
    settings.referenceFilename = GenomicConsensusExperimentalTests::ChimeraFasta;

    const auto names = Workflow::ReferenceNames(settings);
    EXPECT_EQ(4, names.size());
}

TEST(GenomicConsensusExperimentalTest, reference_names_from_filter_windows)
{
    Settings settings;
    settings.filterWindows =
    {
        {"Barcode0--0_Cluster1_Phase1_NumReads297", {300, 600}},
        {"Barcode0--0_Cluster1_Phase1_NumReads297", {2000, 3000}},
        {"Barcode0--0_Cluster0_Phase2_NumReads92", {500, 600}}
    };

    const auto names = Workflow::ReferenceNames(settings);
    EXPECT_EQ(2, names.size());
}

// -----------------------
// Arrow-specific
// -----------------------

// -----------------------
// Plurality-specific
// -----------------------

// -----------------------
// Poa-specific
// -----------------------


// ##########################################################################################
// ##########################################################################################

//TEST(GenomicConsensusExperimentalTests, load_reference_windows_from_fasta)
//{
//    using Input = PacBio::GenomicConsensus::experimental::Input;
//    using Settings = PacBio::GenomicConsensus::experimental::Settings;
//    using GenomicConsensusExperimentalTests::checkInterval;

//    Settings settings;
//    settings.referenceFilename = tests::DataDir + "/chimera_minimal.fasta";

//    const Input input{settings};
//    const auto windows = input.ReferenceWindows();
//    ASSERT_EQ(28, windows.size());

//    SCOPED_TRACE("load_reference_windows_from_fasta");

//    // Barcode0--0_Cluster1_Phase1_NumReads297
//    checkInterval(windows.at(0).interval, 0, 500);
//    checkInterval(windows.at(1).interval, 500, 1000);
//    checkInterval(windows.at(2).interval, 1000, 1500);
//    checkInterval(windows.at(3).interval, 1500, 2000);
//    checkInterval(windows.at(4).interval, 2000, 2500);
//    checkInterval(windows.at(5).interval, 2500, 3000);
//    checkInterval(windows.at(6).interval, 3000, 3152);

//    // Barcode0--0_Cluster1_Phase1_NumReads297
//    checkInterval(windows.at(7).interval, 0, 500);
//    checkInterval(windows.at(8).interval, 500, 1000);
//    checkInterval(windows.at(9).interval, 1000, 1500);
//    checkInterval(windows.at(10).interval, 1500, 2000);
//    checkInterval(windows.at(11).interval, 2000, 2500);
//    checkInterval(windows.at(12).interval, 2500, 3000);
//    checkInterval(windows.at(13).interval, 3000, 3137);

//    // Barcode0--0_Cluster0_Phase2_NumReads92
//    checkInterval(windows.at(14).interval, 0, 500);
//    checkInterval(windows.at(15).interval, 500, 1000);
//    checkInterval(windows.at(16).interval, 1000, 1500);
//    checkInterval(windows.at(17).interval, 1500, 2000);
//    checkInterval(windows.at(18).interval, 2000, 2500);
//    checkInterval(windows.at(19).interval, 2500, 3000);
//    checkInterval(windows.at(20).interval, 3000, 3402);

//    // Barcode0--0_Cluster1_Phase3_NumReads56
//    checkInterval(windows.at(21).interval, 0, 500);
//    checkInterval(windows.at(22).interval, 500, 1000);
//    checkInterval(windows.at(23).interval, 1000, 1500);
//    checkInterval(windows.at(24).interval, 1500, 2000);
//    checkInterval(windows.at(25).interval, 2000, 2500);
//    checkInterval(windows.at(26).interval, 2500, 3000);
//    checkInterval(windows.at(27).interval, 3000, 3151);
//}

//TEST(GenomicConsensusExperimentalTests, enlarged_windows_from_fasta)
//{
//    using Input = PacBio::GenomicConsensus::experimental::Input;
//    using IPoaModel = PacBio::GenomicConsensus::experimental::IPoaModel;
//    using ReferenceWindow = PacBio::GenomicConsensus::experimental::ReferenceWindow;
//    using Settings = PacBio::GenomicConsensus::experimental::Settings;
//    using GenomicConsensusExperimentalTests::checkInterval;

//    Settings settings;
//    settings.referenceFilename = tests::DataDir + "/chimera_minimal.fasta";

//    const Input input{settings};
//    std::vector<ReferenceWindow> eWindows;
//    for (const auto& win : input.ReferenceWindows()) {
//        eWindows.push_back(IPoaModel::EnlargedWindow(win));
//    ASSERT_EQ(28, eWindows.size());

//    SCOPED_TRACE("enlarged_windows_from_fasta");

//    // Barcode0--0_Cluster1_Phase1_NumReads297
//    checkInterval(eWindows.at(0).interval, 0, 505);
//    checkInterval(eWindows.at(1).interval, 495, 1005);
//    checkInterval(eWindows.at(2).interval, 995, 1505);
//    checkInterval(eWindows.at(3).interval, 1495, 2005);
//    checkInterval(eWindows.at(4).interval, 1995, 2505);
//    checkInterval(eWindows.at(5).interval, 2495, 3005);
//    checkInterval(eWindows.at(6).interval, 2995, 3152);

//    // Barcode0--0_Cluster1_Phase1_NumReads297
//    checkInterval(eWindows.at(7).interval, 0, 505);
//    checkInterval(eWindows.at(8).interval, 495, 1005);
//    checkInterval(eWindows.at(9).interval, 995, 1505);
//    checkInterval(eWindows.at(10).interval, 1495, 2005);
//    checkInterval(eWindows.at(11).interval, 1995, 2505);
//    checkInterval(eWindows.at(12).interval, 2495, 3005);
//    checkInterval(eWindows.at(13).interval, 2995, 3137);

//    // Barcode0--0_Cluster0_Phase2_NumReads92
//    checkInterval(eWindows.at(14).interval, 0, 505);
//    checkInterval(eWindows.at(15).interval, 495, 1005);
//    checkInterval(eWindows.at(16).interval, 995, 1505);
//    checkInterval(eWindows.at(17).interval, 1495, 2005);
//    checkInterval(eWindows.at(18).interval, 1995, 2505);
//    checkInterval(eWindows.at(19).interval, 2495, 3005);
//    checkInterval(eWindows.at(20).interval, 2995, 3402);

//    // Barcode0--0_Cluster1_Phase3_NumReads56
//    checkInterval(eWindows.at(21).interval, 0, 505);
//    checkInterval(eWindows.at(22).interval, 495, 1005);
//    checkInterval(eWindows.at(23).interval, 995, 1505);
//    checkInterval(eWindows.at(24).interval, 1495, 2005);
//    checkInterval(eWindows.at(25).interval, 1995, 2505);
//    checkInterval(eWindows.at(26).interval, 2495, 3005);
//    checkInterval(eWindows.at(27).interval, 2995, 3151);
//}

//TEST(GenomicConsensusExperimentalTests, empty_intervals_from_empty_transcript)
//{
//    using Arrow = PacBio::GenomicConsensus::experimental::Arrow;

//    const std::string transcript;
//    const auto intervals = Arrow::TranscriptIntervals(transcript);
//    EXPECT_TRUE(intervals.empty());
//}

//TEST(GenomicConsensusExperimentalTests, single_interval_from_single_op_transcript)
//{
//    using Arrow = PacBio::GenomicConsensus::experimental::Arrow;

//    const std::string transcript = "MMM";
//    const auto intervals = Arrow::TranscriptIntervals(transcript);
//    ASSERT_EQ(1, intervals.size());

//    const auto& interval = intervals.at(0);
//    EXPECT_EQ(0, interval.Left());
//    EXPECT_EQ(3, interval.Right());
//    EXPECT_EQ(3, interval.Length());
//}

//TEST(GenomicConsensusExperimentalTests, intervals_from_transcript)
//{
//    using Arrow = PacBio::GenomicConsensus::experimental::Arrow;

//    const std::string transcript = "MMMRRDDDD";
//    const auto intervals = Arrow::TranscriptIntervals(transcript);
//    ASSERT_EQ(3, intervals.size());

//    auto interval = intervals.at(0);
//    EXPECT_EQ(0, interval.Left());
//    EXPECT_EQ(3, interval.Right());
//    EXPECT_EQ(3, interval.Length());

//    interval = intervals.at(1);
//    EXPECT_EQ(3, interval.Left());
//    EXPECT_EQ(5, interval.Right());
//    EXPECT_EQ(2, interval.Length());

//    interval = intervals.at(2);
//    EXPECT_EQ(5, interval.Left());
//    EXPECT_EQ(9, interval.Right());
//    EXPECT_EQ(4, interval.Length());
//}

//TEST(GenomicConsensusExperimentalTests, calculate_median_from_odd_size)
//{
//    std::vector<size_t> v{5, 6, 4, 3, 2, 6, 7, 9, 3};
//    // { 2, 3, 3, 4, 5, 6, 6, 7 , 9 }
//    const auto median = PacBio::GenomicConsensus::experimental::Arrow::Median(v);
//    EXPECT_EQ(5, median);  // middle
//}

//TEST(GenomicConsensusExperimentalTests, calculate_median_from_even_size)
//{
//    std::vector<size_t> v{5, 6, 4, 3, 2, 6, 7, 3};
//    // { 2, 3, 3, 4, 5, 6, 6, 7 }
//    const auto median = PacBio::GenomicConsensus::experimental::Arrow::Median(v);
//    EXPECT_EQ(4, median);  // avg of middle two (rounds down)
//}

// clang-format on
