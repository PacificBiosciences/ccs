// Author: Lance Hepler

#include "../UnanimityInternalConfig.h"

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
namespace P6C4NoCov {
namespace {

static constexpr const size_t CONTEXT_NUMBER = 1;
static constexpr const size_t OUTCOME_NUMBER = 2;

static constexpr const double kEps = 0.00505052456472967;
static constexpr const double kInvEps = 1.0 - kEps;

class P6C4NoCov_Model : public ModelConfig
{
public:
    static std::set<std::string> Chemistries() { return {"P6-C4"}; }
    static ModelForm Form() { return ModelForm::SNR; }
    P6C4NoCov_Model(const SNR& snr);
    std::unique_ptr<AbstractRecursor> CreateRecursor(const MappedRead& mr,
                                                     double scoreDiff) const override;
    std::vector<TemplatePosition> Populate(const std::string& tpl) const override;
    std::pair<Data::Read, std::vector<MoveType>> SimulateRead(
        std::default_random_engine* const rng, const std::string& tpl,
        const std::string& readname) const override;

    double ExpectedLLForEmission(MoveType move, const AlleleRep& prev, const AlleleRep& curr,
                                 MomentType moment) const override;

private:
    SNR snr_;
    double ctxTrans_[4][2][4];
};

// TODO(lhepler) comments regarding the CRTP
class P6C4NoCovRecursor : public Recursor<P6C4NoCovRecursor>
{
public:
    P6C4NoCovRecursor(const MappedRead& mr, double scoreDiff, double counterWeight);
    static inline std::vector<uint8_t> EncodeRead(const MappedRead& read);
    inline double EmissionPr(MoveType move, uint8_t emission, const AlleleRep& prev,
                             const AlleleRep& curr) const;
    double UndoCounterWeights(size_t nEmissions) const override;

private:
    double counterWeight_;
    double nLgCounterWeight_;
};

static constexpr const double emissionPmf[3][CONTEXT_NUMBER][OUTCOME_NUMBER] = {
    // 0 (match), 1 (mismatch)
    {{kInvEps, kEps / 3.0}},  // MATCH
    {{1.0, 0.0}},             // BRANCH
    {{0.0, 1.0 / 3.0}}        // STICK
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

static constexpr const double snrRanges[2][4] = {
    {0, 0, 0, 0},     // minimum
    {20, 19, 20, 20}  // maximum
};
// For P6-C4 we cap SNR at 20.0 (19.0 for C); as the training set only went that
// high; extrapolation beyond this cap goes haywire because of the higher-order
// terms in the regression model.  See bug 31423.
P6C4NoCov_Model::P6C4NoCov_Model(const SNR& snr)
    : snr_(ClampSNR(snr, SNR{snrRanges[0]}, SNR{snrRanges[1]}))
{
    for (size_t bp = 0; bp < 4; ++bp) {
        for (const bool hp : {false, true}) {
            const auto& params = P6C4NoCovParams[bp][hp];
            const double snr1 = snr_[bp], snr2 = snr1 * snr1, snr3 = snr2 * snr1;
            double sum = 1.0;

            // re-index ctxTrans_ into MBSD from BDMS order
            //                         0123      0123
            //                                   1302 <-- mapping
            const uint8_t mapping[4] = {1, 3, 0, 2};
            ctxTrans_[bp][hp][mapping[0]] = 1.0;
            for (size_t j = 0; j < 3; ++j) {
                double xb =
                    params[j][0] + params[j][1] * snr1 + params[j][2] * snr2 + params[j][3] * snr3;
                xb = std::exp(xb);
                ctxTrans_[bp][hp][mapping[j + 1]] = xb;
                sum += xb;
            }

            for (size_t j = 0; j < 4; ++j)
                ctxTrans_[bp][hp][j] /= sum;
        }
    }
}

std::vector<TemplatePosition> P6C4NoCov_Model::Populate(const std::string& tpl) const
{
    auto rowFetcher = [this](const NCBI2na prev, const NCBI2na curr) -> const double(&)[4]
    {
        const bool hp = (prev.Data() == curr.Data());  // NA -> 0, AA -> 1
        const double(&params)[4] = ctxTrans_[curr.Data()][hp];
        return params;
    };
    return AbstractPopulater(tpl, rowFetcher);
}

std::unique_ptr<AbstractRecursor> P6C4NoCov_Model::CreateRecursor(const MappedRead& mr,
                                                                  double scoreDiff) const
{
    const double counterWeight = CounterWeight(
        [this](size_t ctx, MoveType m) {
            return ctxTrans_[ctx >> 1][ctx & 1][static_cast<uint8_t>(m)];
        },
        [](size_t, MoveType m) {
            switch (m) {
                case MoveType::MATCH:
                    return kInvEps * std::log(kInvEps) + kEps * std::log(kEps / 3.0);
                case MoveType::STICK:
                    return -std::log(3.0);
                default:
                    break;
            }
            return 0.0;
        },
        8);

    return std::make_unique<P6C4NoCovRecursor>(mr, scoreDiff, counterWeight);
}

double P6C4NoCov_Model::ExpectedLLForEmission(const MoveType move, const AlleleRep& prev,
                                              const AlleleRep& curr, const MomentType moment) const
{
    auto cachedEmissionVisitor = [](const MoveType move, const NCBI2na prev, const NCBI2na curr,
                                    const MomentType moment) -> double {
        const double lgThird = -std::log(3.0);
        if (move == MoveType::MATCH) {
            static constexpr const double probMatch = kInvEps;
            static constexpr const double probMismatch = kEps;
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
    };
    return AbstractExpectedLLForEmission(move, prev, curr, moment, cachedEmissionVisitor);
}

P6C4NoCovRecursor::P6C4NoCovRecursor(const MappedRead& mr, double scoreDiff, double counterWeight)
    : Recursor<P6C4NoCovRecursor>(mr, scoreDiff)
    , counterWeight_{counterWeight}
    , nLgCounterWeight_{-std::log(counterWeight_)}
{
}

std::vector<uint8_t> P6C4NoCovRecursor::EncodeRead(const MappedRead& read)
{
    std::vector<uint8_t> result;
    result.reserve(read.Length());

    for (const char bp : read.Seq) {
        result.emplace_back(EncodeBase(bp));
    }

    return result;
}

double P6C4NoCovRecursor::EmissionPr(const MoveType move, const uint8_t emission,
                                     const AlleleRep& prev, const AlleleRep& curr) const
{
    return AbstractEmissionPr(emissionPmf, move, emission, prev, curr) * counterWeight_;
}

double P6C4NoCovRecursor::UndoCounterWeights(const size_t nEmissions) const
{
    return nLgCounterWeight_ * nEmissions;
}

inline std::pair<Data::SNR, std::vector<TemplatePosition>> P6C4NoCov_InitialiseModel(
    std::default_random_engine* const rng, const std::string& tpl)
{
    Data::SNR snrs{0, 0, 0, 0};
    for (uint8_t i = 0; i < 4; ++i) {
        snrs[i] = std::uniform_real_distribution<double>{snrRanges[0][i], snrRanges[1][i]}(*rng);
    }

    const P6C4NoCov_Model model{snrs};
    std::vector<TemplatePosition> transModel = model.Populate(tpl);

    return {snrs, transModel};
}

BaseData P6C4NoCov_GenerateReadData(std::default_random_engine* const rng, const MoveType state,
                                    const AlleleRep& prev, const AlleleRep& curr)
{
    // distribution is arbitrary at the moment, as
    // PW and IPD are not a covariates of the consensus HMM
    std::uniform_int_distribution<uint8_t> pwDistrib{1, 3};
    std::uniform_int_distribution<uint8_t> ipdDistrib{1, 5};

    std::array<double, 4> baseDist;
    for (size_t i = 0; i < 4; ++i) {
        baseDist[i] = AbstractEmissionPr(emissionPmf, state, i, prev, curr);
    }

    std::discrete_distribution<uint8_t> baseDistrib(baseDist.cbegin(), baseDist.cend());

    const char newBase = Data::detail::NCBI2naToASCIIImpl(baseDistrib(*rng));
    const uint8_t newPw = pwDistrib(*rng);
    const uint8_t newIpd = ipdDistrib(*rng);

    return {newBase, newPw, newIpd};
}

std::pair<Data::Read, std::vector<MoveType>> P6C4NoCov_Model::SimulateRead(
    std::default_random_engine* const rng, const std::string& tpl,
    const std::string& readname) const
{
    return SimulateReadImpl(rng, tpl, readname, P6C4NoCov_InitialiseModel,
                            P6C4NoCov_GenerateReadData);
}

}  // namespace anonymous
}  // namespace P6C4NoCov

REGISTER_MODEL_IMPL(P6C4NoCov)

}  // namespace Consensus
}  // namespace PacBio
