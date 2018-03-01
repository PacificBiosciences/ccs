// Author: Armin Töpfer

#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <pbbam/Accuracy.h>
#include <pbbam/LocalContextFlags.h>

#include <pbcopper/cli/CLI.h>

#include <pacbio/ccs/Consensus.h>
#include <pacbio/ccs/ConsensusSettings.h>
#include <pacbio/data/ReadId.h>
#include <pacbio/data/SubreadResultCounter.h>

using namespace PacBio::CCS;
using Subread = ReadType<ReadId>;

TEST(ConsensusTest, TestReadFilter)
{
    std::vector<Subread> data;
    std::string seq =
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    auto movieName = std::make_shared<std::string>("fakeName");
    LocalContextFlags flags = LocalContextFlags::ADAPTER_BEFORE | LocalContextFlags::ADAPTER_AFTER;
    for (int i = 0; i < 10; i++) {
        data.emplace_back(Subread{
            ReadId(movieName, 1, Interval(0, seq.size())), seq, std::vector<uint8_t>(seq.size(), 0),
            std::vector<uint8_t>(seq.size(), 0), flags, .99, SNR(8, 8, 8, 8), "P6-C4"});
    }

    ConsensusSettings settings(PacBio::CLI::Results(ConsensusSettings::CreateCLI("", "")));
    settings.MinLength = 10;
    settings.MinReadScore = 0.0;

    // Nothing filtered
    SubreadResultCounter counter{};
    auto result = FilterReads(data, settings, &counter);
    EXPECT_EQ(0, counter.FilteredBySize);
    // reset
    counter.Success = 0;

    // All removed
    settings.MinLength = 1000;
    auto result2 = FilterReads(data, settings, &counter);
    EXPECT_EQ(10, counter.FilteredBySize);
    EXPECT_EQ(0, counter.Success);
    counter.FilteredBySize = 0;
    counter.Success = 0;

    // Just one
    auto longSeq = seq + seq + seq;
    settings.MinLength = 10;
    data.emplace_back(Subread{ReadId(movieName, 2, Interval(0, longSeq.size())), longSeq,
                              std::vector<uint8_t>(longSeq.size(), 0),
                              std::vector<uint8_t>(longSeq.size(), 0), flags, .99, SNR(8, 8, 8, 8),
                              "P6-C4"});
    auto result3 = FilterReads(data, settings, &counter);
    EXPECT_EQ(1, counter.FilteredBySize);
}
