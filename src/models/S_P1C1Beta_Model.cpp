// Copyright (c) 2011-2016, Pacific Biosciences of California, Inc.
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

#include <cassert>
#include <cmath>
#include <memory>
#include <random>
#include <stdexcept>

#include <pacbio/consensus/ModelConfig.h>
#include <pacbio/data/Read.h>

#include "../ModelFactory.h"
#include "../Recursor.h"
#include "../Simulator.h"
#include "CounterWeight.h"
#include "HelperFunctions.h"

using namespace PacBio::Data;

namespace PacBio {
namespace Consensus {
namespace S_P1C1Beta {
namespace {

class S_P1C1Beta_Model : public ModelConfig
{
    REGISTER_MODEL(S_P1C1Beta_Model);

public:
    static std::set<std::string> Chemistries() { return {"S/P1-C1/beta"}; }
    static ModelForm Form() { return ModelForm::MARGINAL; }
    S_P1C1Beta_Model(const SNR& snr);
    std::unique_ptr<AbstractRecursor> CreateRecursor(const MappedRead& mr,
                                                     double scoreDiff) const override;
    std::vector<TemplatePosition> Populate(const std::string& tpl) const override;
    std::pair<Data::Read, std::vector<MoveType>> SimulateRead(
        std::default_random_engine* const rng, const std::string& tpl,
        const std::string& readname) const override;
    double ExpectedLLForEmission(MoveType move, uint8_t prev, uint8_t curr,
                                 MomentType moment) const override;

private:
    SNR snr_;
};

REGISTER_MODEL_IMPL(S_P1C1Beta_Model);

class S_P1C1Beta_Recursor : public Recursor<S_P1C1Beta_Recursor>
{
public:
    S_P1C1Beta_Recursor(const MappedRead& mr, double scoreDiff, double counterWeight);
    static inline std::vector<uint8_t> EncodeRead(const MappedRead& read);
    inline double EmissionPr(MoveType move, uint8_t emission, uint8_t prev, uint8_t curr) const;
    virtual double UndoCounterWeights(size_t nEmissions) const;

private:
    double counterWeight_;
    double nLgCounterWeight_;
};

static constexpr const double snrRanges[2][4] = {
    {4.0, 4.0, 4.0, 4.0},         // minimum
    {10.65, 10.65, 10.65, 10.65}  // maximum
};

static constexpr const double emissionPmf[3][8][4] = {
    {
        // matchPmf
        {0.980417570, 0.011537479, 0.005804964, 0.002239987},  // AA
        {0.026122324, 0.972937583, 0.000367796, 0.000572296},  // CC
        {0.002544283, 0.002239375, 0.962042375, 0.033173967},  // GG
        {0.000509814, 0.001489097, 0.094228328, 0.903772761},  // TT
        {0.979840156, 0.012582917, 0.005185205, 0.002391722},  // NA
        {0.015528755, 0.984439781, 7.91000E-07, 3.07000E-05},  // NC
        {0.002667013, 0.002095727, 0.961571053, 0.033666207},  // NG
        {0.000506358, 0.001057035, 0.116124340, 0.882312267}   // NT
    },
    {
        // branchPmf
        {1, 0, 0, 0},  // AA
        {0, 1, 0, 0},  // CC
        {0, 0, 1, 0},  // GG
        {0, 0, 0, 1},  // TT
        {1, 0, 0, 0},  // NA
        {0, 1, 0, 0},  // NC
        {0, 0, 1, 0},  // NG
        {0, 0, 0, 1}   // NT
    },
    {
        // stickPmf
        {0.000000000, 0.254503401, 0.574809968, 0.170686631},  // AA
        {0.399446202, 0.000000000, 0.510664061, 0.089889737},  // CC
        {0.505214805, 0.188597323, 0.000000000, 0.306187872},  // GG
        {0.361855644, 0.132870306, 0.505274050, 0.000000000},  // TT
        {0.000000000, 0.210676350, 0.615161689, 0.174161960},  // NA
        {0.357451562, 0.000000000, 0.473482915, 0.169065523},  // NC
        {0.577147745, 0.169785817, 0.000000000, 0.253066438},  // NG
        {0.446834358, 0.144605809, 0.408559833, 0.000000000}   // NT
    }};

static constexpr const double transProbs[8][4] = {
    // Match, Branch, Stick, Delete
    {0.888913751, 0.021169653, 0.034937054, 0.054979542},  // AA
    {0.835822697, 0.036126801, 0.091992041, 0.036058461},  // CC
    {0.886427657, 0.022596867, 0.039619893, 0.051355584},  // GG
    {0.821252207, 0.072798639, 0.068161389, 0.037787765},  // TT
    {0.857630366, 0.072058988, 0.036435296, 0.033875351},  // NA
    {0.846000625, 0.032981179, 0.076759732, 0.044258463},  // NC
    {0.881462348, 0.042444137, 0.039293952, 0.036799562},  // NG
    {0.879087800, 0.022178294, 0.057073518, 0.041660389}   // NT
};

S_P1C1Beta_Model::S_P1C1Beta_Model(const SNR& snr)
    : snr_(ClampSNR(snr, SNR{snrRanges[0]}, SNR{snrRanges[1]}))
{
}

std::vector<TemplatePosition> S_P1C1Beta_Model::Populate(const std::string& tpl) const
{
    std::vector<TemplatePosition> result;

    if (tpl.empty()) return result;

    uint8_t prev = detail::TranslationTable[static_cast<uint8_t>(tpl[0])];
    if (prev > 3) throw std::invalid_argument("invalid character in sequence!");

    for (size_t i = 1; i < tpl.size(); ++i) {
        const uint8_t curr = detail::TranslationTable[static_cast<uint8_t>(tpl[i])];
        if (curr > 3) throw std::invalid_argument("invalid character in sequence!");
        const auto hpAdd = tpl[i - 1] == tpl[i] ? 0 : 4;
        const auto params = transProbs[curr + hpAdd];
        result.emplace_back(TemplatePosition{
            tpl[i - 1], prev,
            params[0],  // match
            params[1],  // branch
            params[2],  // stick
            params[3]   // deletion
        });
        prev = curr;
    }

    result.emplace_back(TemplatePosition{tpl.back(), prev, 1.0, 0.0, 0.0, 0.0});

    return result;
}

std::unique_ptr<AbstractRecursor> S_P1C1Beta_Model::CreateRecursor(const MappedRead& mr,
                                                                   double scoreDiff) const
{
    const double counterWeight = CounterWeight(
        [](size_t ctx, MoveType m) { return transProbs[ctx][static_cast<uint8_t>(m)]; },
        [](size_t ctx, MoveType m) {
            double r = 0.0;
            for (size_t o = 0; o < 4; ++o) {
                const double p = emissionPmf[static_cast<uint8_t>(m)][ctx][o];
                if (p > 0.0) r += p * std::log(p);
            }
            return r;
        },
        8);

    return std::unique_ptr<AbstractRecursor>(new S_P1C1Beta_Recursor(mr, scoreDiff, counterWeight));
}

double S_P1C1Beta_Model::ExpectedLLForEmission(const MoveType move, const uint8_t prev,
                                               const uint8_t curr, const MomentType moment) const
{
    const uint8_t hpAdd = prev == curr ? 0 : 4;
    const uint8_t row = curr + hpAdd;
    double expectedLL = 0;
    for (size_t i = 0; i < 4; i++) {
        double curProb = emissionPmf[static_cast<uint8_t>(move)][row][i];
        double lgCurProb = std::log(curProb);
        if (moment == MomentType::FIRST)
            expectedLL += curProb * lgCurProb;
        else if (moment == MomentType::SECOND)
            expectedLL += curProb * (lgCurProb * lgCurProb);
    }
    return expectedLL;
}

S_P1C1Beta_Recursor::S_P1C1Beta_Recursor(const MappedRead& mr, double scoreDiff,
                                         double counterWeight)
    : Recursor<S_P1C1Beta_Recursor>(mr, scoreDiff)
    , counterWeight_{counterWeight}
    , nLgCounterWeight_{-std::log(counterWeight_)}
{
}

std::vector<uint8_t> S_P1C1Beta_Recursor::EncodeRead(const MappedRead& read)
{
    std::vector<uint8_t> result;
    result.reserve(read.Length());

    for (const char bp : read.Seq) {
        result.emplace_back(EncodeBase(bp));
    }

    return result;
}

inline double S_P1C1Beta_EmissionPr(const MoveType move, const uint8_t emission, const uint8_t prev,
                                    const uint8_t curr)
{
    assert(move != MoveType::DELETION);
    const uint8_t hpAdd = (prev == curr ? 0 : 4);
    // Which row do we want?
    const uint8_t row = curr + hpAdd;
    return emissionPmf[static_cast<uint8_t>(move)][row][emission];
}

double S_P1C1Beta_Recursor::EmissionPr(const MoveType move, const uint8_t emission,
                                       const uint8_t prev, const uint8_t curr) const
{
    return S_P1C1Beta_EmissionPr(move, emission, prev, curr) * counterWeight_;
}

double S_P1C1Beta_Recursor::UndoCounterWeights(const size_t nEmissions) const
{
    return nLgCounterWeight_ * nEmissions;
}

inline std::pair<Data::SNR, std::vector<TemplatePosition>> S_P1C1Beta_InitialiseModel(
    std::default_random_engine* const rng, const std::string& tpl)
{
    Data::SNR snrs{0, 0, 0, 0};
    for (uint8_t i = 0; i < 4; ++i) {
        snrs[i] = std::uniform_real_distribution<double>{snrRanges[0][i], snrRanges[1][i]}(*rng);
    }

    const S_P1C1Beta_Model model{snrs};
    std::vector<TemplatePosition> transModel = model.Populate(tpl);

    return {snrs, transModel};
}

BaseData S_P1C1Beta_GenerateReadData(std::default_random_engine* const rng, const MoveType state,
                                     const uint8_t prev, const uint8_t curr)
{
    static constexpr const std::array<char, 4> bases{{'A', 'C', 'G', 'T'}};

    // distribution is arbitrary at the moment, as
    // PW and IPD are not a covariates of the consensus HMM
    std::uniform_int_distribution<uint8_t> pwDistrib{1, 3};
    std::uniform_int_distribution<uint8_t> ipdDistrib{1, 5};

    std::array<double, 4> baseDist;
    for (size_t i = 0; i < 4; ++i) {
        baseDist[i] = S_P1C1Beta_EmissionPr(state, i, prev, curr);
    }

    std::discrete_distribution<uint8_t> baseDistrib(baseDist.cbegin(), baseDist.cend());

    const char newBase = bases[baseDistrib(*rng)];
    const uint8_t newPw = pwDistrib(*rng);
    const uint8_t newIpd = ipdDistrib(*rng);

    return {newBase, newPw, newIpd};
}

std::pair<Data::Read, std::vector<MoveType>> S_P1C1Beta_Model::SimulateRead(
    std::default_random_engine* const rng, const std::string& tpl,
    const std::string& readname) const
{
    return SimulateReadImpl(rng, tpl, readname, S_P1C1Beta_InitialiseModel,
                            S_P1C1Beta_GenerateReadData);
}

}  // namespace anonymous
}  // namespace S_P1C1Beta
}  // namespace Consensus
}  // namespace PacBio
