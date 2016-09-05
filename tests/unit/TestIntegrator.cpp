// Copyright (c) 2011-2014, Pacific Biosciences of California, Inc.
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

// Author: David Alexander

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <iostream>
#include <limits>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include <pacbio/consensus/MonoMolecularIntegrator.h>
#include <pacbio/consensus/MultiMolecularIntegrator.h>
#include <pacbio/consensus/Mutation.h>
#include <pacbio/data/Sequence.h>

#include "Mutations.h"
#include "RandomDNA.h"

using std::string;
using std::vector;
using std::cerr;
using std::cout;
using std::endl;

using namespace PacBio::Consensus;  // NOLINT
using namespace PacBio::Data;       // NOLINT

using ::testing::UnorderedElementsAreArray;

namespace {

#if EXTENSIVE_TESTING
constexpr int numSamples = 333;
#else
constexpr int numSamples = 3;
#endif

const double prec = 0.001;  // alpha/beta mismatch tolerance
const SNR snr(10, 7, 5, 11);
const string P6C4 = "P6-C4";
const string SP1C1 = "S/P1-C1.1";
const string SP1C1v2 = "S/P1-C1.2";

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
const IntegratorConfig cfg(std::numeric_limits<double>::quiet_NaN());

Read MkRead(const string& seq, const SNR& snr, const string& mdl, const vector<uint8_t>& pw)
{
    vector<uint8_t> ipd(seq.length(), 0);
    return Read("NA", seq, ipd, pw, snr, mdl);
}

#if EXTENSIVE_TESTING
TEST(IntegratorTest, TestLongTemplate)
{
    // TODO: Write a test for a longer molecule
    const string mdl = P6C4;
    vector<uint8_t> pw(longRead.length(), 1);
    MonoMolecularIntegrator ai(longTpl, cfg, snr, mdl);
    EXPECT_EQ(State::VALID,
              ai.AddRead(MappedRead(MkRead(longRead, snr, mdl, pw), StrandType::FORWARD, 0,
                                    longTpl.length(), true, true)));
    EXPECT_NEAR(-148.92614949338801011, ai.LL(), prec);
}

void TestTiming(const string& mdl)
{
    const vector<uint8_t> pws(longRead.length(), 2);
    const size_t nsamp = 5000;
    MonoMolecularIntegrator ai(longTpl, cfg, snr, mdl);
    const auto stime = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < nsamp; ++i)
        EXPECT_EQ(State::VALID,
                  ai.AddRead(MappedRead(MkRead(longRead, snr, mdl, pws), StrandType::FORWARD, 0,
                                        longTpl.length(), true, true)));
    const auto etime = std::chrono::high_resolution_clock::now();
    const auto duration =
        std::chrono::duration_cast<std::chrono::microseconds>(etime - stime).count();
    // std::cout << "avg duration: " << duration / nsamp << "us" << std::endl;
    EXPECT_LT(duration / nsamp, 1500);
}

// disable this test under debug builds (which are not fast enough to pass these timings)
#ifndef NDEBUG
TEST(IntegratorTest, DISABLED_TestLongTemplateTimingP6C4)
#else
TEST(IntegratorTest, TestLongTemplateTimingP6C4)
#endif
{
    TestTiming(P6C4);
}

// disable this test under debug builds (which are not fast enough to pass these timings)
#ifndef NDEBUG
TEST(IntegratorTest, DISABLED_TestLongTemplateTimingSP1C1)
#else
TEST(IntegratorTest, TestLongTemplateTimingSP1C1)
#endif
{
    TestTiming(SP1C1);
}

// disable this test under debug builds (which are not fast enough to pass these timings)
#ifndef NDEBUG
TEST(IntegratorTest, DISABLED_TestLongTemplateTimingSP1C1v2)
#else
TEST(IntegratorTest, TestLongTemplateTimingSP1C1v2)
#endif
{
    TestTiming(SP1C1v2);
}
#endif

std::tuple<string, StrandType> Mutate(const string& tpl, const size_t nmut, std::mt19937* const gen)
{
    string result;

    if (nmut == 0)
        result = tpl;
    else {
        vector<Mutation> muts;
        std::uniform_int_distribution<size_t> rand(0, tpl.length() - 1);
        std::set<size_t> sites;

        while (sites.size() < nmut)
            sites.insert(rand(*gen));

        for (const size_t site : sites) {
            vector<Mutation> possible = Mutations(tpl, site, site + 1);
            std::uniform_int_distribution<size_t> rand2(0, possible.size() - 1);
            muts.push_back(possible[rand2(*gen)]);
        }

        result = ApplyMutations(tpl, &muts);
    }

    std::bernoulli_distribution coin(0.5);

    if (coin(*gen)) return std::make_tuple(result, StrandType::FORWARD);

    return std::make_tuple(ReverseComplement(result), StrandType::REVERSE);
}

template <typename F, typename G>
void MutationEquivalence(const size_t nsamp, const size_t nmut, const F& makeIntegrator,
                         const G& addRead, const string& mdl)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    // increase the floor by 3 because we do not support templates or reads with lt 2 bases,
    //   and we explore edit-1 space around the template to generate reads
    std::uniform_int_distribution<size_t> rand(3 + nmut, 30);

    // count how bad we do
    size_t ntests = 0;
    size_t nerror = 0;

    for (size_t i = 0; i < nsamp; ++i) {
        const string tpl = RandomDNA(rand(gen), &gen);
        vector<Mutation> mutations = Mutations(tpl);
        for (const auto& mut : mutations) {
            string read;
            StrandType strand;
            vector<Mutation> muts{mut};
            const string app = ApplyMutations(tpl, &muts);  // template with mutation applied
            std::tie(read, strand) = Mutate(app, nmut, &gen);
            vector<uint8_t> pws = RandomPW(read.length(), &gen);
            try {
                auto ai1 = makeIntegrator(tpl);
                const auto arr1 = addRead(ai1, MappedRead(MkRead(read, snr, mdl, pws), strand, 0,
                                                          tpl.length(), true, true));
                EXPECT_EQ(State::VALID, arr1);
                if (arr1 != State::VALID) {
                    std::cerr << std::endl << "!! alpha/beta mismatch:" << std::endl;
                    std::cerr << "  " << tpl.length() << ", " << tpl << std::endl;
                    std::cerr << "  " << read.length() << ", " << read << std::endl;
                }
                auto ai2 = makeIntegrator(app);
                const auto arr2 = addRead(ai2, MappedRead(MkRead(read, snr, mdl, pws), strand, 0,
                                                          app.length(), true, true));
                EXPECT_EQ(State::VALID, arr2);
                if (arr2 != State::VALID) {
                    std::cerr << std::endl << "!! alpha/beta mismatch:" << std::endl;
                    std::cerr << "  " << app.length() << ", " << app << std::endl;
                    std::cerr << "  " << read.length() << ", " << read << std::endl;
                }
                const double exp = ai2.LL();
                const double obs0 = ai1.LL();
                const double obs1 = ai1.LL(mut);
                EXPECT_EQ(string(ai1), tpl);
                ai1.ApplyMutations(&muts);
                const double obs2 = ai1.LL();
                // Nigel says the following is not necessarily true
                // if (nmut == 0) EXPECT_LT(obs0, exp);
                // EXPECT_NEAR(obs1, exp, prec);
                // EXPECT_NEAR(obs2, exp, prec);
                EXPECT_EQ(string(ai1), app);
                EXPECT_EQ(string(ai2), app);
                const double diff1 = std::abs(1.0 - obs1 / exp);
                const double diff2 = std::abs(1.0 - obs2 / exp);
                if (diff1 >= prec || diff2 >= prec) {
                    std::cerr << std::endl
                              << "!! intolerable difference: exp: " << exp << ", obs1: " << obs1
                              << ", obs2: " << obs2 << std::endl;
                    std::cerr << "  " << mut << std::endl;
                    std::cerr << "  " << tpl.length() << ", " << tpl << std::endl;
                    std::cerr << "  " << app.length() << ", " << app << std::endl;
                    std::cerr << "  " << ai1.TemplateLength() << ", " << string(ai1) << std::endl;
                    std::stringstream result;
                    std::copy(pws.begin(), pws.end(), std::ostream_iterator<int>(result, " "));
                    std::cerr << "  " << read.length() << ", " << read << " - " << result.str()
                              << std::endl;

                    ++nerror;
                }
            } catch (const std::exception& e) {
                std::cerr << std::endl
                          << "!! caught unexpected exception: " << e.what() << std::endl;
                std::cerr << "  " << mut << std::endl;
                std::cerr << "  " << tpl.length() << ", " << tpl << std::endl;
                std::cerr << "  " << app.length() << ", " << app << std::endl;
                std::cerr << "  " << read.length() << ", " << read << std::endl;
                ++nerror;
            }
            ++ntests;
        }
    }

    EXPECT_EQ(nerror, 0);
    // EXPECT_LT(nerror, ntests / 1000);
}

void MonoEquivalence(const string& mdl)
{
    auto makeMono = [mdl](const string& tpl) {
        return MonoMolecularIntegrator(tpl, cfg, snr, mdl);
    };
    auto monoRead = [](MonoMolecularIntegrator& ai, const MappedRead& mr) {
        return ai.AddRead(mr);
    };
    MutationEquivalence(numSamples, 2, makeMono, monoRead, mdl);
    MutationEquivalence(numSamples, 1, makeMono, monoRead, mdl);
    MutationEquivalence(numSamples, 0, makeMono, monoRead, mdl);
}

TEST(IntegratorTest, TestMonoMutationEquivalenceP6C4) { MonoEquivalence(P6C4); }
TEST(IntegratorTest, TestMonoMutationEquivalenceSP1C1) { MonoEquivalence(SP1C1); }
TEST(IntegratorTest, TestMonoMutationEquivalenceSP1C1v2) { MonoEquivalence(SP1C1v2); }
void MultiEquivalence(const string& mdl)
{
    auto makeMulti = [](const string& tpl) { return MultiMolecularIntegrator(tpl, cfg); };
    auto multiRead = [](MultiMolecularIntegrator& ai, const MappedRead& mr) {
        return ai.AddRead(mr);
    };
    MutationEquivalence(numSamples, 2, makeMulti, multiRead, mdl);
    MutationEquivalence(numSamples, 1, makeMulti, multiRead, mdl);
    MutationEquivalence(numSamples, 0, makeMulti, multiRead, mdl);
}

TEST(IntegratorTest, TestMultiMutationEquivalenceP6C4) { MultiEquivalence(P6C4); }
TEST(IntegratorTest, TestMultiMutationEquivalenceSP1C1) { MultiEquivalence(SP1C1); }
TEST(IntegratorTest, TestMultiMutationEquivalenceSP1C1v2) { MultiEquivalence(SP1C1v2); }
// TODO(lhepler): test multi/mono equivalence
// TODO(lhepler): test multiple mutation testing mono and multi

TEST(IntegratorTest, TestP6C4NoCovAgainstCSharpModel)
{
    const string tpl = "ACGTCGT";
    auto mdl = P6C4;
    MultiMolecularIntegrator ai(tpl, cfg);

    const string readSeq = "ACGTACGT";
    const vector<uint8_t> pw(readSeq.length(), 1);
    EXPECT_EQ(State::VALID,
              ai.AddRead(MappedRead(MkRead("ACGTACGT", snr, mdl, pw), StrandType::FORWARD, 0,
                                    tpl.length(), true, true)));
    auto score = [&ai](Mutation&& mut) { return ai.LL(mut) - ai.LL(); };
    EXPECT_NEAR(-4.74517984808494, ai.LL(), prec);
    EXPECT_NEAR(4.00250386364592, score(Mutation(MutationType::INSERTION, 4, 'A')), prec);
    EXPECT_NEAR(-5.19526526492876, score(Mutation(MutationType::SUBSTITUTION, 2, 'C')), prec);
    EXPECT_NEAR(-4.33430539094949, score(Mutation(MutationType::DELETION, 4)), prec);
    EXPECT_NEAR(-9.70299447206563, score(Mutation(MutationType::DELETION, 6)), prec);
    EXPECT_NEAR(-10.5597017942167, score(Mutation(MutationType::DELETION, 0)), prec);
    EXPECT_NEAR(-0.166992912601578, score(Mutation(MutationType::SUBSTITUTION, 4, 'A')), prec);
    EXPECT_NEAR(-1.60697112438296, score(Mutation(MutationType::INSERTION, 4, 'G')), prec);
}

TEST(IntegratorTest, TestFailAddRead)
{
    const string tpl = "A";
    const vector<uint8_t> pw(tpl.length(), 1);
    auto mdl = P6C4;
    MultiMolecularIntegrator ai(tpl, cfg);

    EXPECT_EQ(State::TEMPLATE_TOO_SMALL,
              ai.AddRead(MappedRead(MkRead(tpl, snr, mdl, pw), StrandType::FORWARD, 0, tpl.length(),
                                    true, true)));
}

TEST(IntegratorTest, TestSuccessAddRead)
{
    const string tpl = "AA";
    const vector<uint8_t> pw(tpl.length(), 1);
    const auto mdl = P6C4;
    MultiMolecularIntegrator ai(tpl, cfg);

    EXPECT_EQ(State::VALID, ai.AddRead(MappedRead(MkRead(tpl, snr, mdl, pw), StrandType::FORWARD, 0,
                                                  tpl.length(), true, true)));
}

}  // namespace anonymous
