// Author: Lance Hepler

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iterator>
#include <set>
#include <string>
#include <vector>

using std::string;
using std::vector;

#include <sys/stat.h>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <pacbio/consensus/Integrator.h>
#include <pacbio/consensus/ModelSelection.h>
#include <pacbio/data/Read.h>
#include <pacbio/data/State.h>

#include "TestData.h"

using namespace PacBio::Consensus;  // NOLINT
using namespace PacBio::Data;       // NOLINT

namespace LoadModelsTests {
const SNR snr(10, 7, 5, 11);

const string longTpl =
    "GGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTTCTTCGTCAT"
    "AACTTAATGTTTTTATTTAAAATACCCTCTGAAAAGAAAGGAAACGACAGGTGCTGAAAGCGAGCTTTTTGGCC"
    "TCTGTCGTTTCCTTTCTCTGTTTTTGTCCGTGGAATGAACAATGGAAGTCAACAAAAAGCAGCTGGCTGACATT"
    "TTCGGTGCGAGTATCCGTACCATTCAGAACTGGCAGGAACAGGGAATGCCCGTTCTGCGAGGCGGTGGCAAGGG"
    "TAATGAGGTGCTTTATGACTCTGCCGCCGTCATAAAATGGTATGCCGAAAGGGATGCTGAAATTGAGAACGAAA"
    "AGCTGCGCCGGGAGGTTGAAGAACTGCGGCAGGCCAGCGAGGCAGATCTCCAGCCAGGAACTATTGAGTACGAA"
    "CGCCATCGACTTACGCGTGCGCAGGCCGACGCACAGGAACTGAAGAATGCCAGAGACTCCGCTGAAGTGGTGGA"
    "AACCGCATTCTGTACTTTCGTGCTGTCGCGGATCGCAGGTGAAATTGCCAGTATTCTCGACGGGCTCCCCCTGT"
    "CGGTGCAGCGGCGTTTTCCGGAACTGGAAAACCGACATGTTGATTTCCTGAAACGGGATATCATCAAAGCCATG"
    "AACAAAGCAGCCGCGCTGGATGAACTGATACCGGGGTTGCTGAGTGAATATATCGAACAGTCAGGTTAACAGGC"
    "TGCGGCATTTTGTCCGCGCCGGGCTTCGCTCACTGTTCAGGCCGGAGCCACAGACCGCCGTTGAATGGGCGGAT"
    "GCTAATTACTATCTCCCGAAAGAATC";

const string longRead =
    "GGGCGGCGACCTCGCGGGTTTTCGCTATTTCTGAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTTCTTCGTCAT"
    "AACTTAATGTTTTTATTTAAAATACCCTCTGAAAAGAAAGGAAACGACAGGTGCTGAAAGCGAGCTTTTTGGCC"
    "TCTGTCGTTTCCTTTCTCTGTTTTTGTCCGTGGAATGAACAATGGAAGTCAACAAAAAGCAGCTGGCTGACATT"
    "TTCGGTGGAGTATCCGTACCATTCAGAACTGGCAGGACAGGGAATGCCCGTTCTGCGAGGCGGTGGCAAGGGTA"
    "ATGAGGTGCTTTATGACTCTGCCGCCGTCATAAAATGGTATGCCGAAAGGGATGCTGAAATTGAGAACGAATAG"
    "CTGCGCCGGGAGGTTGAAGAACTGCGGCAGGCCAGCGAGGCAGATCTCCAGCCAGGAACTATTGAGTACGAACG"
    "CCATCGACTTACGCGTGCGCAGGCCGACGCACAGGAACTGAAGAATGCCAGAGACTCCGCTGAAGTGGTGGAAA"
    "CCGCATTCCCCTGTACTTTCGTGCTGTCGCGGATCGCAGGTGAAATTGCCAGTATTCTCGACGGGCTCCCCCTG"
    "TCGGTGCAGCGGCGTTTTCCGGAACTGGAAAACCGACATGTTGATTTCCTGAAACGGGATATCATCAAAGCCAT"
    "GAACAAAGCAGCCGCGCTGGATGAACTGATACCGGGGTTGCTGAGTGAATATATCGAACAGTCAGGTTAACAGG"
    "CTGCGGCATTTTGTCCGCGCCGGGCTTCGCTCACTGTTCAGGCCGGAGCCACAGACCGCCGTTGAACGGATGCT"
    "AATTACTATCTCCCGAAAGAATC";

const std::vector<uint8_t> longPws(longRead.size(), 10);
const IntegratorConfig cfg(-100.0);  // disable zscore filtering

Read MkRead(const string& seq, const SNR& snr, const string& mdl, const vector<uint8_t>& pw)
{
    vector<uint8_t> ipd(seq.length(), 0);
    return Read("NA", seq, ipd, pw, snr, mdl);
}
}  // namespace LoadModelsTests

TEST(LoadModelsTest, SupportedChemistries)
{
    const std::set<std::string> chem = {"P6-C4",     "S/P1-C1/beta", "S/P1-C1.1",  "S/P1-C1.2",
                                        "S/P1-C1.3", "S/P2-C2",      "S/P2-C2/5.0"};
    const std::set<std::string> supp = SupportedChemistries();
    std::set<std::string> diff;
    std::set_difference(chem.begin(), chem.end(), supp.begin(), supp.end(),
                        std::inserter(diff, diff.begin()));
    EXPECT_TRUE(diff.size() == 0);
}

TEST(LoadModelsTest, Malformed)
{
    std::string malformed = tests::DataDir + "/Malformed.json";
    struct stat st;
    ASSERT_TRUE(stat(malformed.c_str(), &st) == 0 && S_ISREG(st.st_mode));
    ASSERT_FALSE(LoadModels(malformed));
}

// disable this test because we cannot load S/P1-C1.2 more than once
// TEST(LoadModelsTest, DISABLED_SingleFile)
// {
//     std::string sp1c1v2 = tests::DataDir + "/arrow/SP1C1v2.json";
//     ASSERT_TRUE(LoadModelFromFile(sp1c1v2));
//     std::set<std::string> chems = SupportedChemistries();
//     ASSERT_TRUE(chems.find("S/P1-C1.2::PwSnr") != chems.end());
// }

TEST(LoadModelsTest, Directory)
{
    ASSERT_TRUE(LoadModels(tests::DataDir + "/arrow"));
    std::set<std::string> chems = SupportedModels();
    ASSERT_TRUE(chems.find("S/P1-C1/beta::Marginal::FromFile") != chems.end());
    ASSERT_TRUE(chems.find("S/P1-C1.1::PwSnrA::FromFile") != chems.end());
    ASSERT_TRUE(chems.find("S/P1-C1.2::PwSnr::FromFile") != chems.end());
    ASSERT_TRUE(chems.find("S/P2-C2/5.0::PwSnr::FromFile") != chems.end());

    // test identity between S/P1-C1/beta and S/P1-C1/beta::Marginal (loaded)
    {
        Integrator ai1(LoadModelsTests::longTpl, LoadModelsTests::cfg);
        EXPECT_EQ(State::VALID,
                  ai1.AddRead(MappedRead(
                      LoadModelsTests::MkRead(LoadModelsTests::longRead, LoadModelsTests::snr,
                                              "S/P1-C1/beta::Marginal::Compiled",
                                              LoadModelsTests::longPws),
                      StrandType::FORWARD, 0, LoadModelsTests::longTpl.length(), true, true)));

        Integrator ai2(LoadModelsTests::longTpl, LoadModelsTests::cfg);
        EXPECT_EQ(State::VALID,
                  ai2.AddRead(MappedRead(
                      LoadModelsTests::MkRead(LoadModelsTests::longRead, LoadModelsTests::snr,
                                              "S/P1-C1/beta::Marginal::FromFile",
                                              LoadModelsTests::longPws),
                      StrandType::FORWARD, 0, LoadModelsTests::longTpl.length(), true, true)));

        EXPECT_NEAR(ai1.LL(), ai2.LL(), 1.0e-5);
    }

    // test identity between S/P1-C1.1 and S/P1-C1.1::PwSnrA (loaded)
    {
        Integrator ai1(LoadModelsTests::longTpl, LoadModelsTests::cfg);
        EXPECT_EQ(
            State::VALID,
            ai1.AddRead(MappedRead(
                LoadModelsTests::MkRead(LoadModelsTests::longRead, LoadModelsTests::snr,
                                        "S/P1-C1.1::PwSnrA::Compiled", LoadModelsTests::longPws),
                StrandType::FORWARD, 0, LoadModelsTests::longTpl.length(), true, true)));

        Integrator ai2(LoadModelsTests::longTpl, LoadModelsTests::cfg);
        EXPECT_EQ(
            State::VALID,
            ai2.AddRead(MappedRead(
                LoadModelsTests::MkRead(LoadModelsTests::longRead, LoadModelsTests::snr,
                                        "S/P1-C1.1::PwSnrA::FromFile", LoadModelsTests::longPws),
                StrandType::FORWARD, 0, LoadModelsTests::longTpl.length(), true, true)));

        EXPECT_NEAR(ai1.LL(), ai2.LL(), 1.0e-5);
    }

    // test identity between S/P1-C1.2 and S/P1-C1.2::PwSnr
    {
        Integrator ai1(LoadModelsTests::longTpl, LoadModelsTests::cfg);
        EXPECT_EQ(
            State::VALID,
            ai1.AddRead(MappedRead(
                LoadModelsTests::MkRead(LoadModelsTests::longRead, LoadModelsTests::snr,
                                        "S/P1-C1.2::PwSnr::Compiled", LoadModelsTests::longPws),
                StrandType::FORWARD, 0, LoadModelsTests::longTpl.length(), true, true)));

        Integrator ai2(LoadModelsTests::longTpl, LoadModelsTests::cfg);
        EXPECT_EQ(
            State::VALID,
            ai2.AddRead(MappedRead(
                LoadModelsTests::MkRead(LoadModelsTests::longRead, LoadModelsTests::snr,
                                        "S/P1-C1.2::PwSnr::FromFile", LoadModelsTests::longPws),
                StrandType::FORWARD, 0, LoadModelsTests::longTpl.length(), true, true)));

        EXPECT_NEAR(ai1.LL(), ai2.LL(), 1.0e-5);
    }

    // test identity between S/P1-C1.2 and S/P1-C1.2::PwSnr
    {
        Integrator ai1(LoadModelsTests::longTpl, LoadModelsTests::cfg);
        EXPECT_EQ(State::VALID,
                  ai1.AddRead(MappedRead(
                      LoadModelsTests::MkRead(LoadModelsTests::longRead, LoadModelsTests::snr,
                                              "S/P2-C2::PwSnr::Compiled", LoadModelsTests::longPws),
                      StrandType::FORWARD, 0, LoadModelsTests::longTpl.length(), true, true)));

        Integrator ai2(LoadModelsTests::longTpl, LoadModelsTests::cfg);
        EXPECT_EQ(
            State::VALID,
            ai2.AddRead(MappedRead(
                LoadModelsTests::MkRead(LoadModelsTests::longRead, LoadModelsTests::snr,
                                        "S/P1-C1.2::PwSnr::FromFile", LoadModelsTests::longPws),
                StrandType::FORWARD, 0, LoadModelsTests::longTpl.length(), true, true)));

        EXPECT_NEAR(ai1.LL(), ai2.LL(), 1.0e-5);
    }

    // test identity between S/P2-C2/5.0 and S/P2-C2/5.0::PwSnr
    {
        Integrator ai1(LoadModelsTests::longTpl, LoadModelsTests::cfg);
        EXPECT_EQ(
            State::VALID,
            ai1.AddRead(MappedRead(
                LoadModelsTests::MkRead(LoadModelsTests::longRead, LoadModelsTests::snr,
                                        "S/P2-C2/5.0::PwSnr::Compiled", LoadModelsTests::longPws),
                StrandType::FORWARD, 0, LoadModelsTests::longTpl.length(), true, true)));

        Integrator ai2(LoadModelsTests::longTpl, LoadModelsTests::cfg);
        EXPECT_EQ(
            State::VALID,
            ai2.AddRead(MappedRead(
                LoadModelsTests::MkRead(LoadModelsTests::longRead, LoadModelsTests::snr,
                                        "S/P2-C2/5.0::PwSnr::FromFile", LoadModelsTests::longPws),
                StrandType::FORWARD, 0, LoadModelsTests::longTpl.length(), true, true)));

        EXPECT_NEAR(ai1.LL(), ai2.LL(), 1.0e-5);
    }
}

TEST(LoadModelsTest, UpdateBundle)
{
    const char* varname = "SMRT_CHEMISTRY_BUNDLE_DIR";
    EXPECT_EQ(0, setenv(varname, tests::DataDir.c_str(), 0));

    std::set<std::string> chems = SupportedModels();
    ASSERT_TRUE(chems.find("S/P1-C1/beta::Marginal::Bundled") != chems.end());
    ASSERT_TRUE(chems.find("S/P1-C1.1::PwSnrA::Bundled") != chems.end());
    ASSERT_TRUE(chems.find("S/P1-C1.2::PwSnr::Bundled") != chems.end());

    EXPECT_EQ(0, unsetenv(varname));
}

#if EXTENSIVE_TESTING
#ifndef NDEBUG
TEST(LoadModelsTest, DISABLED_ModelTiming)
#else
TEST(LoadModelsTest, ModelTiming)
#endif
{
    // load required models just in case they haven't been already
    LoadModels(tests::DataDir + "/arrow");

    const size_t nsamp = 100;
    const std::vector<std::string> mdls = {
        "S/P1-C1/beta::Marginal::FromFile", "S/P1-C1.1::PwSnrA::FromFile",
        "S/P1-C1.2::PwSnr::FromFile", "S/P2-C2/5.0::PwSnr::FromFile"};
    for (const auto mdl : mdls) {
        Integrator ai(LoadModelsTests::longTpl, LoadModelsTests::cfg);
        const auto stime = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < nsamp; ++i)
            EXPECT_EQ(State::VALID,
                      ai.AddRead(MappedRead(
                          LoadModelsTests::MkRead(LoadModelsTests::longRead, LoadModelsTests::snr,
                                                  mdl, LoadModelsTests::longPws),
                          StrandType::FORWARD, 0, LoadModelsTests::longTpl.length(), true, true)));
        const auto etime = std::chrono::high_resolution_clock::now();
        const auto duration =
            std::chrono::duration_cast<std::chrono::microseconds>(etime - stime).count();
        // std::cout << mdl << " avg duration: " << duration / nsamp << "us" << std::endl;
        EXPECT_LT(duration / nsamp, 1500);
    }
}
#endif

TEST(LoadModelsTest, ModelOverride)
{
    Integrator ai1(LoadModelsTests::longTpl, LoadModelsTests::cfg);
    EXPECT_EQ(State::VALID,
              ai1.AddRead(MappedRead(
                  LoadModelsTests::MkRead(LoadModelsTests::longRead, LoadModelsTests::snr,
                                          "S/P1-C1.2", LoadModelsTests::longPws),
                  StrandType::FORWARD, 0, LoadModelsTests::longTpl.length(), true, true)));

    ASSERT_TRUE(OverrideModel("S/P1-C1.2"));

    Integrator ai2(LoadModelsTests::longTpl, LoadModelsTests::cfg);
    EXPECT_EQ(State::VALID,
              ai2.AddRead(MappedRead(
                  LoadModelsTests::MkRead(LoadModelsTests::longRead, LoadModelsTests::snr,
                                          "S/P1-C1.1", LoadModelsTests::longPws),
                  StrandType::FORWARD, 0, LoadModelsTests::longTpl.length(), true, true)));

    EXPECT_NEAR(ai1.LL(), ai2.LL(), 1.0e-5);

    ASSERT_TRUE(UnOverrideModel());
}
