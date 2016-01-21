
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

double SP1C1NoCovParams[8][4] =
    { // Match, Branch, Stick, Delete
    {0.869389092, 0.027611193, 0.027662784, 0.075336931}, // AA
    {0.838429016, 0.035022031, 0.08333472 , 0.043214233}, // CC
    {0.864109073, 0.030418552, 0.036400241, 0.069072134}, // GG
    {0.767901703, 0.069269032, 0.100084066, 0.062745199}, // TT
    {0.848057757, 0.073970928, 0.039412757, 0.038558558}, // NA
    {0.904383595, 0.019235302, 0.067664074, 0.00871703},  // NC
    {0.844087036, 0.072763801, 0.041130542, 0.04201862},  // NG
    {0.846995852, 0.012811015, 0.075218145, 0.064974988}, // NT
    };

SP1C1NoCovModel::SP1C1NoCovModel(const SNR& snr) : snr_(snr), counterWeight_{1.0}
{

}

std::vector<TemplatePosition> SP1C1NoCovModel::Populate(const std::string& tpl) const
{
    std::vector<TemplatePosition> result;

    if (tpl.empty()) return result;

    for (size_t i = 1; i < tpl.size(); ++i) {
        const uint8_t b = detail::TranslationTable[static_cast<uint8_t>(tpl[i])];

        if (b > 3) throw std::invalid_argument("invalid character in sequence!");

        const bool hpAdd = tpl[i - 1] == tpl[i] ? 4 : 0;
        const auto params = SP1C1NoCovParams[b + hpAdd];
        
        result.emplace_back(TemplatePosition{
            tpl[i - 1],
            params[0],  // match
            params[1],  // branch
            params[2],  // stick
            params[3]   // deletion
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

    constexpr double pr = 0.04091752;

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
