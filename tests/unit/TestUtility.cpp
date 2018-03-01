// Authors: David Alexander, Lance Hepler

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <seqan/seq_io.h>

#include <pacbio/io/Utility.h>

#include "TestUtility.h"

bool LoadFastaSequences(std::string fastaFname, std::vector<std::string>& ids,
                        std::vector<std::string>& seqs)
{
    seqan::StringSet<seqan::CharString> sIds;
    seqan::StringSet<seqan::CharString> sSeqs;

    seqan::SeqFileIn seqFileIn(fastaFname.c_str());
    readRecords(sIds, sSeqs, seqFileIn);

    for (auto& sId : sIds) {
        ids.emplace_back(toCString(sId));
    }
    for (auto& sSeq : sSeqs) {
        seqs.emplace_back(toCString(sSeq));
    }

    return true;
}