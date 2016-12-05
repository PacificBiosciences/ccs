// Copyright (c) 2016, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: Lance Hepler

#include <algorithm>
#include <chrono>
#include <iterator>
#include <set>
#include <string>
#include <vector>

#include <sys/stat.h>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <pacbio/consensus/ModelSelection.h>
#include <pacbio/consensus/MonoMolecularIntegrator.h>
#include <pacbio/data/Read.h>
#include <pacbio/data/State.h>

#include "TestData.h"

using namespace std;
using namespace PacBio::Consensus;  // NOLINT
using namespace PacBio::Data;       // NOLINT

namespace {
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

const std::vector<uint8_t> longPws(longRead.size(), 2);
const IntegratorConfig cfg(-100.0);  // disable zscore filtering

Read MkRead(const string& seq, const SNR& snr, const string& mdl, const vector<uint8_t>& pw)
{
    vector<uint8_t> ipd(seq.length(), 0);
    return Read("NA", seq, ipd, pw, snr, mdl);
}
}

TEST(LoadModelsTest, SupportedChemistries)
{
    const std::set<std::string> chem = {"P6-C4",     "S/P1-C1/beta", "S/P1-C1.1",
                                        "S/P1-C1.2", "S/P1-C1.3",    "S/P2-C2"};
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
//     std::string sp1c1v2 = tests::DataDir + "/params/SP1C1v2.json";
//     ASSERT_TRUE(LoadModelFromFile(sp1c1v2));
//     std::set<std::string> chems = SupportedChemistries();
//     ASSERT_TRUE(chems.find("S/P1-C1.2::PwSnr") != chems.end());
// }

TEST(LoadModelsTest, Directory)
{
    ASSERT_TRUE(LoadModels(tests::DataDir + "/params"));
    std::set<std::string> chems = SupportedModels();
    ASSERT_TRUE(chems.find("S/P1-C1/beta::Marginal::FromFile") != chems.end());
    ASSERT_TRUE(chems.find("S/P1-C1.1::PwSnrA::FromFile") != chems.end());
    ASSERT_TRUE(chems.find("S/P1-C1.2::PwSnr::FromFile") != chems.end());
    ASSERT_TRUE(chems.find("S/P2-C2::PwSnr::FromFile") != chems.end());

// test identity between S/P1-C1/beta and S/P1-C1/beta::Marginal (loaded)
//   disabled until S_P1C1Beta is fixed
#if 0
    {
        MonoMolecularIntegrator ai1(longTpl, cfg, snr, "S/P1-C1/beta::Compiled");
        EXPECT_EQ(State::VALID,
                  ai1.AddRead(MappedRead(MkRead(longRead, snr, "S/P1-C1/beta::Compiled", longPws),
                                         StrandType::FORWARD, 0, longTpl.length(), true, true)));

        MonoMolecularIntegrator ai2(longTpl, cfg, snr, "S/P1-C1/beta::Marginal");
        EXPECT_EQ(State::VALID,
                  ai2.AddRead(MappedRead(MkRead(longRead, snr, "S/P1-C1/beta::Marginal", longPws),
                                         StrandType::FORWARD, 0, longTpl.length(), true, true)));

        EXPECT_NEAR(ai1.LL(), ai2.LL(), 1.0e-5);
    }
#endif

    // test identity between S/P1-C1.1 and S/P1-C1.1::PwSnrA (loaded)
    {
        MonoMolecularIntegrator ai1(longTpl, cfg, snr, "S/P1-C1.1::PwSnrA::Compiled");
        EXPECT_EQ(State::VALID, ai1.AddRead(MappedRead(
                                    MkRead(longRead, snr, "S/P1-C1.1::PwSnrA::Compiled", longPws),
                                    StrandType::FORWARD, 0, longTpl.length(), true, true)));

        MonoMolecularIntegrator ai2(longTpl, cfg, snr, "S/P1-C1.1::PwSnrA::FromFile");
        EXPECT_EQ(State::VALID, ai2.AddRead(MappedRead(
                                    MkRead(longRead, snr, "S/P1-C1.1::PwSnrA::FromFile", longPws),
                                    StrandType::FORWARD, 0, longTpl.length(), true, true)));

        EXPECT_NEAR(ai1.LL(), ai2.LL(), 1.0e-5);
    }

    // test identity between S/P1-C1.2 and S/P1-C1.2::PwSnr
    {
        MonoMolecularIntegrator ai1(longTpl, cfg, snr, "S/P1-C1.2::PwSnr::Compiled");
        EXPECT_EQ(State::VALID, ai1.AddRead(MappedRead(
                                    MkRead(longRead, snr, "S/P1-C1.2::PwSnr::Compiled", longPws),
                                    StrandType::FORWARD, 0, longTpl.length(), true, true)));

        MonoMolecularIntegrator ai2(longTpl, cfg, snr, "S/P1-C1.2::PwSnr::FromFile");
        EXPECT_EQ(State::VALID, ai2.AddRead(MappedRead(
                                    MkRead(longRead, snr, "S/P1-C1.2::PwSnr::FromFile", longPws),
                                    StrandType::FORWARD, 0, longTpl.length(), true, true)));

        EXPECT_NEAR(ai1.LL(), ai2.LL(), 1.0e-5);
    }

    // test identity between S/P1-C1.2 and S/P1-C1.2::PwSnr
    {
        MonoMolecularIntegrator ai1(longTpl, cfg, snr, "S/P2-C2::PwSnr::Compiled");
        EXPECT_EQ(State::VALID,
                  ai1.AddRead(MappedRead(MkRead(longRead, snr, "S/P2-C2::PwSnr::Compiled", longPws),
                                         StrandType::FORWARD, 0, longTpl.length(), true, true)));

        MonoMolecularIntegrator ai2(longTpl, cfg, snr, "S/P2-C2::PwSnr::FromFile");
        EXPECT_EQ(State::VALID,
                  ai2.AddRead(MappedRead(MkRead(longRead, snr, "S/P2-C2::PwSnr::FromFile", longPws),
                                         StrandType::FORWARD, 0, longTpl.length(), true, true)));

        EXPECT_NEAR(ai1.LL(), ai2.LL(), 1.0e-5);
    }
}

#if EXTENSIVE_TESTING
#ifndef NDEBUG
TEST(LoadModelsTest, DISABLED_ModelTiming)
#else
TEST(LoadModelsTest, ModelTiming)
#endif
{
    // load required models just in case they haven't been already
    LoadModels(tests::DataDir + "/params");

    const size_t nsamp = 100;
    const std::vector<std::string> mdls = {
        "S/P1-C1/beta::Marginal::FromFile", "S/P1-C1.1::PwSnrA::FromFile",
        "S/P1-C1.2::PwSnr::FromFile", "S/P2-C2::PwSnr::FromFile"};
    for (const auto mdl : mdls) {
        MonoMolecularIntegrator ai(longTpl, cfg, snr, mdl);
        const auto stime = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < nsamp; ++i)
            EXPECT_EQ(State::VALID,
                      ai.AddRead(MappedRead(MkRead(longRead, snr, mdl, longPws),
                                            StrandType::FORWARD, 0, longTpl.length(), true, true)));
        const auto etime = std::chrono::high_resolution_clock::now();
        const auto duration =
            std::chrono::duration_cast<std::chrono::microseconds>(etime - stime).count();
        // std::cout << "avg duration: " << duration / nsamp << "us" << std::endl;
        EXPECT_LT(duration / nsamp, 1500);
    }
}
#endif

TEST(LoadModelsTest, ModelOverride)
{
    MonoMolecularIntegrator ai1(longTpl, cfg, snr, "S/P1-C1.2");
    EXPECT_EQ(State::VALID,
              ai1.AddRead(MappedRead(MkRead(longRead, snr, "S/P1-C1.2", longPws),
                                     StrandType::FORWARD, 0, longTpl.length(), true, true)));

    ASSERT_TRUE(OverrideModel("S/P1-C1.2"));

    MonoMolecularIntegrator ai2(longTpl, cfg, snr, "S/P1-C1.1");
    EXPECT_EQ(State::VALID,
              ai2.AddRead(MappedRead(MkRead(longRead, snr, "S/P1-C1.1", longPws),
                                     StrandType::FORWARD, 0, longTpl.length(), true, true)));

    EXPECT_NEAR(ai1.LL(), ai2.LL(), 1.0e-5);

    ASSERT_TRUE(UnOverrideModel());
}
