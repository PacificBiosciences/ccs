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
#include <stdexcept>

#include <pacbio/consensus/ModelConfig.h>
#include <pacbio/consensus/Read.h>

#include "../ModelFactory.h"
#include "../Recursor.h"

namespace PacBio {
namespace Consensus {
namespace {

class SP1C1BetaNoCovModel : public ModelConfig
{
    REGISTER_MODEL(SP1C1BetaNoCovModel);

public:
    static std::set<std::string> Names() { return {"S/P1-C1/beta"}; }
    SP1C1BetaNoCovModel(const SNR& snr);
    std::unique_ptr<AbstractRecursor> CreateRecursor(std::unique_ptr<AbstractTemplate>&& tpl,
                                                     const MappedRead& mr, double scoreDiff) const;
    std::vector<TemplatePosition> Populate(const std::string& tpl) const;
    double ExpectedLogLikelihoodForMatchEmission(uint8_t prev, uint8_t curr, bool secondMoment) const;
    double ExpectedLogLikelihoodForStickEmission(uint8_t prev, uint8_t curr, bool secondMoment) const;
    double ExpectedLogLikelihoodForBranchEmission(uint8_t prev, uint8_t curr, bool secondMoment) const;

private:
    SNR snr_;
};

REGISTER_MODEL_IMPL(SP1C1BetaNoCovModel);

class SP1C1BetaNoCovRecursor : public Recursor<SP1C1BetaNoCovRecursor>
{
public:
    SP1C1BetaNoCovRecursor(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr,
                           double scoreDiff);
    static inline std::vector<uint8_t> EncodeRead(const MappedRead& read);
    static inline double EmissionPr(MoveType move, uint8_t emission, uint8_t prev, uint8_t curr);
    virtual double UndoCounterWeights(size_t nEmissions) const;
};

double emissionPmf[3][8][4] = {{
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

double SP1C1BetaNoCovParams[8][4] = {
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

SP1C1BetaNoCovModel::SP1C1BetaNoCovModel(const SNR& snr) : snr_(snr) {}
std::vector<TemplatePosition> SP1C1BetaNoCovModel::Populate(const std::string& tpl) const
{
    std::vector<TemplatePosition> result;

    if (tpl.empty()) return result;

    uint8_t prev = detail::TranslationTable[static_cast<uint8_t>(tpl[0])];
    if (prev > 3) throw std::invalid_argument("invalid character in sequence!");

    for (size_t i = 1; i < tpl.size(); ++i) {
        const uint8_t curr = detail::TranslationTable[static_cast<uint8_t>(tpl[i])];
        if (curr > 3) throw std::invalid_argument("invalid character in sequence!");
        const bool hpAdd = tpl[i - 1] == tpl[i] ? 0 : 4;
        const auto params = SP1C1BetaNoCovParams[curr + hpAdd];
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

std::unique_ptr<AbstractRecursor> SP1C1BetaNoCovModel::CreateRecursor(
    std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr, double scoreDiff) const
{
    return std::unique_ptr<AbstractRecursor>(new SP1C1BetaNoCovRecursor(
        std::forward<std::unique_ptr<AbstractTemplate>>(tpl), mr, scoreDiff));
}

inline double CalculateExpectedLogLikelihoodOfOutcomeRow(const int index, const uint8_t prev, const uint8_t curr, const bool secondMoment)  {
        const uint8_t hpAdd = prev == curr ? 0 : 4;
        const uint8_t row = curr + hpAdd;
        double expectedLL = 0;
        for(size_t i = 0; i < 4; i++) {
            double curProb = emissionPmf[index][row][i];
            double lgCurProb = std::log(curProb);
            if(!secondMoment) {
                expectedLL +=  curProb * lgCurProb;
            } else {
                expectedLL += curProb * pow(lgCurProb, 2.0);
            }
        }
        return expectedLL;
    }
    

double SP1C1BetaNoCovModel::ExpectedLogLikelihoodForMatchEmission(uint8_t prev, uint8_t curr, bool secondMoment) const {
        return CalculateExpectedLogLikelihoodOfOutcomeRow(static_cast<uint8_t>(MoveType::MATCH), prev, curr, secondMoment);
}
double SP1C1BetaNoCovModel::ExpectedLogLikelihoodForStickEmission(uint8_t prev, uint8_t curr, bool secondMoment) const {
        return CalculateExpectedLogLikelihoodOfOutcomeRow(static_cast<uint8_t>(MoveType::STICK), prev, curr, secondMoment);
}
double SP1C1BetaNoCovModel::ExpectedLogLikelihoodForBranchEmission(uint8_t prev, uint8_t curr, bool secondMoment) const {
        return CalculateExpectedLogLikelihoodOfOutcomeRow(static_cast<uint8_t>(MoveType::BRANCH), prev, curr, secondMoment);
}

SP1C1BetaNoCovRecursor::SP1C1BetaNoCovRecursor(std::unique_ptr<AbstractTemplate>&& tpl,
                                               const MappedRead& mr, double scoreDiff)
    : Recursor<SP1C1BetaNoCovRecursor>(std::forward<std::unique_ptr<AbstractTemplate>>(tpl), mr,
                                       scoreDiff)
{
}

std::vector<uint8_t> SP1C1BetaNoCovRecursor::EncodeRead(const MappedRead& read)
{
    std::vector<uint8_t> result;

    for (const char bp : read.Seq) {
        const uint8_t em = detail::TranslationTable[static_cast<uint8_t>(bp)];
        if (em > 3) throw std::invalid_argument("invalid character in read!");
        result.emplace_back(em);
    }

    return result;
}

double SP1C1BetaNoCovRecursor::EmissionPr(MoveType move, const uint8_t emission, const uint8_t prev,
                                          const uint8_t curr)
{
    assert(move != MoveType::DELETION);
    const uint8_t hpAdd = prev == curr ? 0 : 4;
    // Which row do we want?
    const uint8_t row = curr + hpAdd;
    return emissionPmf[static_cast<uint8_t>(move)][row][emission];
}

double SP1C1BetaNoCovRecursor::UndoCounterWeights(const size_t nEmissions) const { return 0; }
}  // namespace anonymous
}  // namespace Consensus
}  // namespace PacBio
