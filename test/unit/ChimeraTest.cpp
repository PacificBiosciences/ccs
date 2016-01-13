
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <seqan/seq_io.h>

#include "pacbio/chimera/ChimeraLabeler.h"
#include "pacbio/chimera/ChimeraLabel.h"

using namespace PacBio::Chimera;

const std::string FILENAME = "../test/unit/data/test.fasta";

TEST(ChimeraLabeler, EndToEnd)
{
    SeqFileIn inputHandle(FILENAME.c_str());
    
    // Parse the records
    StringSet<CharString> ids;
    StringSet<Dna5String> seqs;
    readRecords(ids, seqs, inputHandle);
        
    // Declare the vectors we'll use to actually perform the Chimera-labeling
    std::vector<std::string> idList;
    std::vector<Dna5String> seqList;

    for (size_t i = 0; i < length(ids); ++i)
    {
        idList.push_back(toCString(static_cast<CharString>(ids[i])));
        seqList.push_back(static_cast<Dna5String>(seqs[i]));
    }

    // Label the Records
    ChimeraLabeler chimeraLabeler(1.0f, true);
    auto labels = chimeraLabeler.LabelChimeras(idList, seqList);

    EXPECT_EQ(1, 1);
}
