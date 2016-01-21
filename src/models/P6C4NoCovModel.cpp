
#include <cassert>
#include <cmath>
#include <stdexcept>

#include <pacbio/consensus/ModelConfig.h>

#include "../ModelFactory.h"

namespace PacBio {
namespace Consensus {
namespace {

class P6C4NoCovModel : public ModelConfig
{
    REGISTER_MODEL(P6C4NoCovModel);

public:
    P6C4NoCovModel(const SNR& snr);
    std::vector<TemplatePosition> Populate(const std::string& tpl) const;
    double BaseEmissionPr(MoveType move, char from, char to) const;
    double CovEmissionPr(MoveType move, uint8_t cov, const char from, const char to) const;
    double UndoCounterWeights(size_t nEmissions) const;

    static std::string Name() { return "P6-C4"; }
private:
    SNR snr_;
    double counterWeight_;
};

REGISTER_MODEL_IMPL(P6C4NoCovModel);

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

P6C4NoCovModel::P6C4NoCovModel(const SNR& snr) : snr_(snr), counterWeight_{1.0}
{
    constexpr auto bases = "ACGT";
    constexpr MoveType moves[] = {MoveType::MATCH, MoveType::BRANCH, MoveType::STICK};

    double baseEmission = 0.0;
    double covEmission = 0.0;

    for (size_t m = 0; m < 3; ++m) {
        for (size_t i = 0; i < 4; ++i)
            for (size_t j = 0; j < 4; ++j)
                baseEmission += BaseEmissionPr(moves[m], bases[i], bases[j]);

        for (uint8_t c = 0; c < 20; ++c)
            covEmission += CovEmissionPr(moves[m], c, 0, 0);
    }

    baseEmission /= (3 * 4 * 4);
    covEmission /= (3 * 20);
    counterWeight_ /= (baseEmission * covEmission);
}

std::vector<TemplatePosition> P6C4NoCovModel::Populate(const std::string& tpl) const
{
    std::vector<TemplatePosition> result;

    if (tpl.empty()) return result;

    for (size_t i = 1; i < tpl.size(); ++i) {
        const uint8_t b = detail::TranslationTable[static_cast<uint8_t>(tpl[i])];

        if (b > 3) throw std::invalid_argument("invalid character in sequence!");

        const bool hp = tpl[i - 1] == tpl[i];  // NA -> 0, AA -> 1
        const auto params = P6C4NoCovParams[b][hp];
        const double snr = snr_[b], snr2 = snr * snr, snr3 = snr2 * snr;
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
            tpl[i - 1],
            tprobs[1],  // match
            1.0 / sum,  // branch
            tprobs[2],  // stick
            tprobs[0]   // deletion
        });
    }

    result.emplace_back(TemplatePosition{tpl.back(), 0.0, 0.0, 0.0, 0.0});

    return result;
}

double P6C4NoCovModel::BaseEmissionPr(MoveType move, const char from, const char to) const
{
    assert(move != MoveType::DELETION);

    if (move == MoveType::BRANCH)
        return 1.0;
    else if (move == MoveType::STICK)
        return 1.0 / 3.0;

    constexpr double pr = 0.00505052456472967;

    if (from == to) return 1.0 - pr;

    return pr / 3.0;
}

double P6C4NoCovModel::CovEmissionPr(MoveType move, const uint8_t, const char from, const char to) const
{
    assert(move != MoveType::DELETION);
    return 1.0 * counterWeight_;
}

double P6C4NoCovModel::UndoCounterWeights(const size_t nEmissions) const
{
    return -std::log(counterWeight_) * nEmissions;
}

}  // namespace anonymous
}  // namespace Consensus
}  // namespace PacBio
