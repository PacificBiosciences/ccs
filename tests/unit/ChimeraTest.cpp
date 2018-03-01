// Author: Armin Töpfer

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <seqan/seq_io.h>

#include "pacbio/chimera/ChimeraLabel.h"
#include "pacbio/chimera/ChimeraLabeler.h"

#include "TestData.h"

using namespace PacBio::Chimera;

const std::string FILENAME1 = tests::DataDir + "/chimera_minimal.fasta";
const std::string FILENAME2 = tests::DataDir + "/chimera_extensive.fasta";

TEST(ChimeraLabeler, MinimalEndToEnd)
{
    using namespace seqan;

    SeqFileIn inputHandle(FILENAME1.c_str());

    // Parse the records
    StringSet<CharString> ids;
    StringSet<Dna5String> seqs;
    readRecords(ids, seqs, inputHandle);

    // Declare the vectors we'll use to actually perform the Chimera-labeling
    std::vector<std::string> idList;
    std::vector<std::string> seqList;

    for (size_t i = 0; i < length(ids); ++i) {
        idList.emplace_back(toCString(static_cast<CharString>(ids[i])));
        seqList.emplace_back(toCString(static_cast<CharString>(seqs[i])));
    }

    // Label the Records
    ChimeraLabeler chimeraLabeler;
    const auto& labels = chimeraLabeler.LabelChimeras(idList, seqList);

    // Expect 4 labels
    EXPECT_EQ(labels.size(), 4);

    // Expect the first 2 (non-Chimeric) sequences to have no score
    EXPECT_EQ(labels[0].score, -1.0);
    EXPECT_EQ(labels[1].score, -1.0);

    // Expect the next (homologous) sequence to have a low, non-zero score
    EXPECT_LT(labels[2].score, 1.0);
    EXPECT_GT(labels[2].score, 0.0);

    // Expect the last (Chimeric) sequence to have a high score
    EXPECT_GT(labels[3].score, 1.0);
}

#if EXTENSIVE_TESTING
TEST(ChimeraLabeler, ExtensiveEndToEnd)
{
    using namespace seqan;

    SeqFileIn inputHandle(FILENAME2.c_str());

    // Parse the records
    StringSet<CharString> ids;
    StringSet<Dna5String> seqs;
    readRecords(ids, seqs, inputHandle);

    // Declare the vectors we'll use to actually perform the Chimera-labeling
    std::vector<std::string> idList;
    std::vector<std::string> seqList;

    for (size_t i = 0; i < length(ids); ++i) {
        idList.push_back(toCString(static_cast<CharString>(ids[i])));
        seqList.push_back(toCString(static_cast<CharString>(seqs[i])));
    }

    // Label the Records
    ChimeraLabeler chimeraLabeler;
    const auto& labels = chimeraLabeler.LabelChimeras(idList, seqList);

    // Expect 10 labels
    EXPECT_EQ(labels.size(), 10);

    // Expect the first 6 (non-Chimeric) sequences to have low scores
    EXPECT_LT(labels[0].score, 1.0);
    EXPECT_LT(labels[1].score, 1.0);
    EXPECT_LT(labels[2].score, 1.0);
    EXPECT_LT(labels[3].score, 1.0);
    EXPECT_LT(labels[4].score, 1.0);
    EXPECT_LT(labels[5].score, 1.0);

    // Expect the last 4 (Chimeric) sequences to have high scores
    EXPECT_GT(labels[6].score, 1.0);
    EXPECT_GT(labels[7].score, 1.0);
    EXPECT_GT(labels[8].score, 1.0);
    EXPECT_GT(labels[9].score, 1.0);
}
#endif
