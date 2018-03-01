// Author: Lance Hepler

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <iostream>
#include <string>
#include <tuple>

using std::string;

#include <pacbio/consensus/Integrator.h>
#include <pacbio/consensus/Polish.h>
#include <pacbio/data/Read.h>
#include <pacbio/data/Sequence.h>

using namespace PacBio::Consensus;  // NOLINT
using namespace PacBio::Data;       // NOLINT

namespace PolishTests {

Read MkRead(const std::string& seq, const SNR& snr, const std::string& mdl)
{
    std::vector<uint8_t> cov(seq.length(), 0);
    return Read("NA", seq, cov, cov, snr, mdl);
}

const SNR snr(10, 7, 5, 11);
const string mdl("P6-C4");

TEST(PolishTest, Basic)
{
    Integrator ai("GCGTCGT", IntegratorConfig());

    ai.AddRead(MappedRead(MkRead("ACGTACGT", snr, mdl), StrandType::FORWARD, 0, 7, true, true));
    ai.AddRead(MappedRead(MkRead(ReverseComplement("ACGACGT"), snr, mdl), StrandType::REVERSE, 0, 7,
                          true, true));
    ai.AddRead(MappedRead(MkRead("ACGACGT", snr, mdl), StrandType::FORWARD, 0, 7, true, true));

    const auto result = Polish(&ai, PolishConfig());

    EXPECT_TRUE(result.hasConverged);
    EXPECT_EQ("ACGACGT", string(ai));
}

TEST(PolishTest, DiTriRepeat)
{
    //                       1  2  31 2 3
    const string tpl = "ACGTCAGCAGCAGAGAGTGCA";
    //                       1  2  3  41 2 3 4
    const string read = "ACGTCAGCAGCAGCAGAGAGAGTGCA";
    Integrator ai(tpl, IntegratorConfig());

    ai.AddRead(
        MappedRead(MkRead(read, snr, mdl), StrandType::FORWARD, 0, tpl.length(), true, true));
    ai.AddRead(MappedRead(MkRead(ReverseComplement(read), snr, mdl), StrandType::REVERSE, 0,
                          tpl.length(), true, true));
    ai.AddRead(
        MappedRead(MkRead(read, snr, mdl), StrandType::FORWARD, 0, tpl.length(), true, true));
    // ai.AddRead(MappedRead(MkRead(ReverseComplement(tpl), snr, mdl), StrandType::REVERSE, 0, tpl.length(), true, true));
    // ai.AddRead(MappedRead(MkRead(tpl, snr, mdl), StrandType::FORWARD, 0, tpl.length(), true, true));

    if (false) {
        const auto noresult = Polish(&ai, PolishConfig());

        EXPECT_TRUE(noresult.hasConverged);
        EXPECT_EQ(tpl, string(ai));
    }

    const auto result = PolishRepeats(&ai, RepeatConfig());

    EXPECT_TRUE(result.hasConverged);
    EXPECT_EQ(read, string(ai));
}
}  // namespace PolishTests
