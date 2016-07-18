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
#include <pacbio/data/Read.h>

#include "../ModelFactory.h"
#include "../Recursor.h"

using namespace PacBio::Data;

namespace PacBio {
namespace Consensus {
namespace {

constexpr double kEps = 0.00505052456472967;
constexpr double kCounterWeight = 1.894736842105264607;

class P6C4NoCovModel : public ModelConfig
{
    REGISTER_MODEL(P6C4NoCovModel);

public:
    static std::set<std::string> Names() { return {"P6-C4"}; }
    P6C4NoCovModel(const SNR& snr);
    std::unique_ptr<AbstractRecursor> CreateRecursor(std::unique_ptr<AbstractTemplate>&& tpl,
                                                     const MappedRead& mr, double scoreDiff) const;
    std::vector<TemplatePosition> Populate(const std::string& tpl) const;
    double ExpectedLLForEmission(MoveType move, uint8_t prev, uint8_t curr,
                                 MomentType moment) const;

private:
    SNR snr_;
};

REGISTER_MODEL_IMPL(P6C4NoCovModel);

// TODO(lhepler) comments regarding the CRTP
class P6C4NoCovRecursor : public Recursor<P6C4NoCovRecursor>
{
public:
    P6C4NoCovRecursor(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr,
                      double scoreDiff);
    static inline std::vector<uint8_t> EncodeRead(const MappedRead& read);
    static inline double EmissionPr(MoveType move, uint8_t emission, uint8_t prev, uint8_t curr);
    virtual double UndoCounterWeights(size_t nEmissions) const;
};

double P6C4NoCovParams[4][2][3][4] = {
    { // A
     {// NA
      {2.35936060895653, -0.463630601682986, 0.0179206897766131, -0.000230839937063052},
      {3.22847830625841, -0.0886820214931539, 0.00555981712798726, -0.000137686231186054},
      {-0.101031042923432, -0.0138783767832632, -0.00153408019582419, 7.66780338484727e-06}},
     {// AA
      {3.76122480667588, -0.536010820176981, 0.0275375059387171, -0.000470200724345621},
      {3.57517725358548, -0.0257545295375707, -0.000163673803286944, 5.3256984681724e-06},
      {0.858421613302247, -0.0276654216841666, -8.85549766507732e-05, -4.85355908595337e-05}}},
    { // C
     {// NC
      {5.956054206161, -1.71886470811695, 0.153315470604752, -0.00474488595513198},
      {3.89418464416296, -0.174182841558867, 0.0171719290275442, -0.000653629721359769},
      {2.40532887070852, -0.652606650098156, 0.0688783864119339, -0.00246479494650594}},
     {// CC
      {5.66725538674764, -1.10462196933913, 0.0879811093908922, -0.00259393800835979},
      {4.11682756767018, -0.124758322644639, 0.00659795177909886, -0.000361914629195461},
      {3.17103818507405, -0.729020290806687, 0.0749784690396837, -0.00262779517495421}}},
    { // G
     {// NG
      {3.53508304630569, -0.788027301381263, 0.0469367803413207, -0.00106221924705805},
      {2.85440184222226, 0.166346531056167, -0.0166161828155307, 0.000439492705370092},
      {0.238188180807376, 0.0589443522886522, -0.0123401045958974, 0.000336854126836293}},
     {// GG
      {3.81920778703052, -0.540309003502589, 0.0389569264893982, -0.000901245733796236},
      {3.31322216145728, 0.123514009118836, -0.00807401406655071, 0.000230843924466035},
      {2.06006877520527, -0.451486652688621, 0.0375212898173045, -0.000937676250926241}}},
    { // T
     {// NT
      {5.36199280681367, -1.46099908985536, 0.126755291030074, -0.0039102734460725},
      {3.41597143103046, -0.066984162951578, 0.0138944877787003, -0.000558939998921912},
      {1.37371376794871, -0.246963827944892, 0.0209674231346363, -0.000684856715039738}},
     {// TT
      {5.39308368236762, -1.32931568057267, 0.107844580241936, -0.00316462903462847},
      {4.21031404956015, -0.347546363361823, 0.0293839179303896, -0.000893802212450644},
      {2.33143889851302, -0.586068444099136, 0.040044954697795, -0.000957298861394191}}}};

// For P6-C4 we cap SNR at 20.0 (19.0 for C); as the training set only went that
// high; extrapolation beyond this cap goes haywire because of the higher-order
// terms in the regression model.  See bug 31423.
P6C4NoCovModel::P6C4NoCovModel(const SNR& snr)
    : snr_(ClampSNR(snr, SNR(0, 0, 0, 0), SNR(20, 19, 20, 20)))
{
}

std::vector<TemplatePosition> P6C4NoCovModel::Populate(const std::string& tpl) const
{
    std::vector<TemplatePosition> result;

    if (tpl.empty()) return result;

    uint8_t prev = detail::TranslationTable[static_cast<uint8_t>(tpl[0])];
    if (prev > 3) throw std::invalid_argument("invalid character in sequence!");

    for (size_t i = 1; i < tpl.size(); ++i) {
        const uint8_t curr = detail::TranslationTable[static_cast<uint8_t>(tpl[i])];
        if (curr > 3) throw std::invalid_argument("invalid character in sequence!");
        const bool hp = tpl[i - 1] == tpl[i];  // NA -> 0, AA -> 1
        const auto params = P6C4NoCovParams[curr][hp];
        const double snr = snr_[curr], snr2 = snr * snr, snr3 = snr2 * snr;
        double tprobs[3];
        double sum = 1.0;

        for (size_t j = 0; j < 3; ++j) {
            double xb =
                params[j][0] + snr * params[j][1] + snr2 * params[j][2] + snr3 * params[j][3];
            xb = std::exp(xb);
            tprobs[j] = xb;
            sum += xb;
        }

        for (size_t j = 0; j < 3; ++j)
            tprobs[j] /= sum;

        result.emplace_back(TemplatePosition{
            tpl[i - 1], prev,
            tprobs[1],  // match
            1.0 / sum,  // branch
            tprobs[2],  // stick
            tprobs[0]   // deletion
        });

        prev = curr;
    }

    result.emplace_back(TemplatePosition{tpl.back(), prev, 1.0, 0.0, 0.0, 0.0});

    return result;
}

std::unique_ptr<AbstractRecursor> P6C4NoCovModel::CreateRecursor(
    std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr, double scoreDiff) const
{
    return std::unique_ptr<AbstractRecursor>(
        new P6C4NoCovRecursor(std::forward<std::unique_ptr<AbstractTemplate>>(tpl), mr, scoreDiff));
}

double P6C4NoCovModel::ExpectedLLForEmission(const MoveType move, const uint8_t prev,
                                             const uint8_t curr, const MomentType moment) const
{
    const double lgThird = -std::log(3.0);
    if (move == MoveType::MATCH) {
        constexpr double probMatch = 1.0 - kEps;
        constexpr double probMismatch = kEps;
        const double lgMatch = std::log(probMatch);
        const double lgMismatch = lgThird + std::log(probMismatch);
        if (moment == MomentType::FIRST)
            return probMatch * lgMatch + probMismatch * lgMismatch;
        else if (moment == MomentType::SECOND)
            return probMatch * (lgMatch * lgMatch) + probMismatch * (lgMismatch * lgMismatch);
    } else if (move == MoveType::BRANCH)
        return 0.0;
    else if (move == MoveType::STICK) {
        if (moment == MomentType::FIRST)
            return lgThird;
        else if (moment == MomentType::SECOND)
            return lgThird * lgThird;
    }
    throw std::invalid_argument("invalid move!");
}

P6C4NoCovRecursor::P6C4NoCovRecursor(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr,
                                     double scoreDiff)
    : Recursor<P6C4NoCovRecursor>(std::forward<std::unique_ptr<AbstractTemplate>>(tpl), mr,
                                  scoreDiff)
{
}

std::vector<uint8_t> P6C4NoCovRecursor::EncodeRead(const MappedRead& read)
{
    std::vector<uint8_t> result;

    for (const char bp : read.Seq) {
        const uint8_t em = detail::TranslationTable[static_cast<uint8_t>(bp)];
        if (em > 3) throw std::invalid_argument("invalid character in read!");
        result.emplace_back(em);
    }

    return result;
}

double P6C4NoCovRecursor::EmissionPr(MoveType move, const uint8_t emission, const uint8_t prev,
                                     const uint8_t curr)
{
    assert(move != MoveType::DELETION);

    // probability of a mismatch
    constexpr double tbl[3][2] = {
        // 0 (match), 1 (mismatch)
        {1.0 - kEps, kEps / 3.0},  // MATCH
        {1.0, 0.0},                // BRANCH
        {0.0, 1.0 / 3.0}           // STICK
    };

    return tbl[static_cast<uint8_t>(move)][curr != emission] * kCounterWeight;
}

double P6C4NoCovRecursor::UndoCounterWeights(const size_t nEmissions) const
{
    return -std::log(kCounterWeight) * nEmissions;
}

}  // namespace anonymous
}  // namespace Consensus
}  // namespace PacBio
