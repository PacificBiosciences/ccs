
#include <cassert>
#include <cmath>
#include <stdexcept>

#include <pacbio/consensus/ModelConfig.h>

#include "../ModelFactory.h"

namespace PacBio {
namespace Consensus {
namespace {

class SP1C1NoCovModel : public ModelConfig
{
    REGISTER_MODEL(SP1C1NoCovModel);

public:
    SP1C1NoCovModel(const SNR& snr);
    std::vector<TemplatePosition> Populate(const std::string& tpl) const;
    double BaseEmissionPr(MoveType move, char from, char to) const;
    double CovEmissionPr(MoveType move, uint8_t cov) const;
    double UndoCounterWeights(size_t nEmissions) const;

    static std::string Name() { return "S/P1-C1"; }
private:
    SNR snr_;
    double counterWeight_;
};

REGISTER_MODEL_IMPL(SP1C1NoCovModel);

double SP1C1NoCovParams[4][2][3][4] = {
    { // A
     {// NA
      {4.65845057808354, -0.415840910444698, -0.11828058244689, 0.00956578673791941},
      {3.57612618844509, 0.468988528255572, -0.137940992879249, 0.0084375333741803},
      {1.72030498433343, -0.492926937734869, -0.0102314103567128, 0.00249347261642999}},
     {// AA
      {-1.05509248823551, 2.84395927296652, -0.538272420181546, 0.0280168485992841},
      {0.443975466103867, 2.15978168706798, -0.366850584349875, 0.0190036316025785},
      {-2.14040849388288, 1.81773844138814, -0.337854936333007, 0.0163975738768652}}},
    { // C
     {// NC
      {3.34486645444278, -0.565986396762248, 0.0228093724792962, -0.00024878685432142},
      {4.93052555254383, -0.381207322234503, 0.0323305306361318, -0.000816152981640072},
      {1.87204897233723, -0.634450686340275, 0.0572254276700295, -0.0014388392076932}},
     {// CC
      {7.91046549090261, -1.75322053199142, 0.152476251650544, -0.00448514638942164},
      {6.54740336306551, -1.01984567386713, 0.105531870279225, -0.00340223372227548},
      {2.21930560154163, -1.09315562577391, 0.125436615376555, -0.00406650091873246}}},
    { // G
     {// NG
      {12.1549366815382, -3.29550224416333, 0.393699153551291, -0.0188299063737343},
      {12.060230629607, -3.15696465161496, 0.474523764795806, -0.0247128009140047},
      {13.5402075889642, -4.91053683674892, 0.68734627463777, -0.0334749656338937}},
     {// GG
      {8.55697621774934, -1.81755140542446, 0.202131631358331, -0.0106290897661297},
      {9.33699464190863, -1.99025580512387, 0.287407757742778, -0.0150179408777114},
      {7.85389951738019, -2.44062223524733, 0.304523048222906, -0.0145051788136893}}},
    { // T
     {// NT
      {-0.801627056990592, 0.352876944067311, -0.0599994635240104, 0.0026873766811099},
      {1.73606606926821, 0.0580083328827675, 0.0175987540988494, -0.000896074886804205},
      {-0.675377454833453, -0.660096907028115, 0.0869371941163534, -0.00284109121410621}},
     {// TT
      {1.76660426955674, -0.0583563240212852, -0.0421223392285456, 0.00234373801558385},
      {2.66316147241342, 0.0713472492667887, -0.0137033054612303, 0.000449808985830026},
      {1.08751469080523, -1.19581914070065, 0.122707565437219, -0.0036860243024887}}}};

SP1C1NoCovModel::SP1C1NoCovModel(const SNR& snr) : snr_(snr), counterWeight_{1.0}
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
            covEmission += CovEmissionPr(moves[m], c);
    }

    baseEmission /= (3 * 4 * 4);
    covEmission /= (3 * 20);
    counterWeight_ /= (baseEmission * covEmission);
}

std::vector<TemplatePosition> SP1C1NoCovModel::Populate(const std::string& tpl) const
{
    std::vector<TemplatePosition> result;

    if (tpl.empty()) return result;

    for (size_t i = 1; i < tpl.size(); ++i) {
        const uint8_t b = detail::TranslationTable[static_cast<uint8_t>(tpl[i])];

        if (b > 3) throw std::invalid_argument("invalid character in sequence!");

        const bool hp = tpl[i - 1] == tpl[i];  // NA -> 0, AA -> 1
        const auto params = SP1C1NoCovParams[b][hp];
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

double SP1C1NoCovModel::BaseEmissionPr(MoveType move, const char from, const char to) const
{
    assert(move != MoveType::DELETION);

    if (move == MoveType::BRANCH)
        return 1.0;
    else if (move == MoveType::STICK)
        return 1.0 / 3.0;

    constexpr double pr = 0.0341694930009249;

    if (from == to) return 1.0 - pr;

    return pr / 3.0;
}

double SP1C1NoCovModel::CovEmissionPr(MoveType move, const uint8_t) const
{
    assert(move != MoveType::DELETION);
    return 1.0 * counterWeight_;
}

double SP1C1NoCovModel::UndoCounterWeights(const size_t nEmissions) const
{
    return -std::log(counterWeight_) * nEmissions;
}

}  // namespace anonymous
}  // namespace Consensus
}  // namespace PacBio
