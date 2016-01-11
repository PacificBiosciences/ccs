
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
      { -0.616100703516245, 0.634945507055802, -0.0989522077243331, 0.00421623115653814  },
      { 2.77631504649577, 0.147901650215569, -0.0154729270843096, 0.00025551668513119  },
      { 0.307586459085395, -0.014990694086659, 0.00749505147412215, -0.000686010705660055  }},
     {// AA
      {-1.05509248823551, 2.84395927296652, -0.538272420181546, 0.0280168485992841},
      {0.443975466103867, 2.15978168706798, -0.366850584349875, 0.0190036316025785},
      {-2.14040849388288, 1.81773844138814, -0.337854936333007, 0.0163975738768652}}},
    { // C
     {// NC
      { 0.388160256039391, -0.0288976239192285, 0.00307318451468241, -0.000234158838210871  },
      { 3.14875514091988, -0.0106823353354246, -0.000401845910921559, -1.25310835565835e-05  },
      { 0.597855536671197, -0.0861810319129838, 0.00740873528590919, -0.000242619026252324  }},
     {// CC
      { -0.303292047753244, 0.152394700979222, -0.0138577900465039, 0.000300102948928903  },
      { 2.71937148455391, 0.0530024828239804, -0.00271845220498571, -4.94530582078834e-06  },
      { 0.472021667180638, -0.0874865552273888, 0.00861464812493361, -0.000217844712414285  } }},
    { // G
     {// NG
      { -0.0165995309100124, 0.123898027084033, -0.0144787041574539, 0.000256376930639249  },
      { 3.57113329561504, -0.167727645988598, 0.0178590733740708, -0.000697471999270054  },
      { 0.264642340709538, -0.0218814584903329, 0.000784672327698186, -2.43997496128624e-05  }},
     {// GG
      { 0.207758700467849, 0.0558652602208208, -0.00964460409277763, 0.000248839532159918  },
      { 2.84517075801878, 0.065554236133839, -0.00752914799222252, 0.000227628713195767  },
      { -1.13392815509634, 0.383932117824245, -0.0413502957513709, 0.00141766336011068  }}},
    { // T
     {// NT
      { -0.929954967132838, 0.2680810142706, -0.0197216063129832, 0.000355976980966589  },
      { 2.31185369504364, 0.163946029245349, -0.0130352658819774, 0.000280578875574133  },
      { -0.540980932507914, 0.172322333326211, -0.0145677914882085, 0.000375186381284634  }},
     {// TT
      { -0.848432413231851, 0.293137416686333, -0.0274167595019264, 0.000699587877966592  },
      { 2.4788614394413, 0.12465919478065, -0.0111013993546001, 0.000269775302903984  },
      { 0.471325777768195, -0.125172007001656, 0.0106653934624698, -0.000269135228378833  }
    }}};

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

    constexpr double pr = 0.107510957040264;

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
