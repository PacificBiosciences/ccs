// Author: Lance Hepler

#include "../UnanimityInternalConfig.h"

#include <cassert>
#include <cmath>
#include <memory>
#include <random>
#include <stdexcept>

#include <pacbio/consensus/ModelConfig.h>
#include <pacbio/data/Read.h>
#include <pacbio/exception/ModelError.h>

#include "../JsonHelpers.h"
#include "../ModelFactory.h"
#include "../ModelFormFactory.h"
#include "../Recursor.h"
#include "../Simulator.h"
#include "CounterWeight.h"
#include "HelperFunctions.h"

using namespace PacBio::Data;

namespace PacBio {
namespace Consensus {
namespace PwSnr {
namespace {

using MalformedModelFile = PacBio::Exception::MalformedModelFile;

static constexpr const size_t CONTEXT_NUMBER = 16;
static constexpr const size_t OUTCOME_NUMBER = 12;

// fwd decl
class PwSnrModelCreator;

class PwSnrModel : public ModelConfig
{
public:
    PwSnrModel(const PwSnrModelCreator* params, const SNR& snr);
    std::unique_ptr<AbstractRecursor> CreateRecursor(const MappedRead& mr,
                                                     double scoreDiff) const override;
    std::vector<TemplatePosition> Populate(const std::string& tpl) const override;
    std::pair<Data::Read, std::vector<MoveType>> SimulateRead(
        std::default_random_engine* const rng, const std::string& tpl,
        const std::string& readname) const override;
    double ExpectedLLForEmission(MoveType move, const AlleleRep& prev, const AlleleRep& curr,
                                 MomentType moment) const override;

    friend class PwSnrInitializeModel;

private:
    double CalculateExpectedLLForEmission(const size_t move, const uint8_t row,
                                          const size_t moment) const;

    const PwSnrModelCreator* params_;
    SNR snr_;
    double ctxTrans_[CONTEXT_NUMBER][4];
    double cachedEmissionExpectations_[CONTEXT_NUMBER][3][2];
};

class PwSnrRecursor : public Recursor<PwSnrRecursor>
{
public:
    PwSnrRecursor(const MappedRead& mr, double scoreDiff, double counterWeight,
                  const PwSnrModelCreator* params);

    static std::vector<uint8_t> EncodeRead(const MappedRead& read);
    double EmissionPr(MoveType move, uint8_t emission, const AlleleRep& prev,
                      const AlleleRep& curr) const;
    double UndoCounterWeights(size_t nEmissions) const override;

private:
    const PwSnrModelCreator* params_;
    double counterWeight_;
    double nLgCounterWeight_;
};

class PwSnrModelCreator : public ModelCreator
{
    friend class PwSnrModel;
    friend class PwSnrRecursor;
    friend class PwSnrInitializeModel;
    friend class PwSnrGenerateReadData;

public:
    static ModelForm Form() { return ModelForm::PWSNR; }
    PwSnrModelCreator(const boost::property_tree::ptree& pt);
    std::unique_ptr<ModelConfig> Create(const SNR& snr) const override
    {
        return std::make_unique<PwSnrModel>(this, snr);
    };

private:
    double snrRanges_[4][2];
    double emissionPmf_[3][CONTEXT_NUMBER][OUTCOME_NUMBER];
    double transitionParams_[CONTEXT_NUMBER][3][4];
};

inline double PwSnrModel::CalculateExpectedLLForEmission(const size_t move, const uint8_t row,
                                                         const size_t moment) const
{
    double expectedLL = 0;
    for (size_t i = 0; i < OUTCOME_NUMBER; i++) {
        double curProb = params_->emissionPmf_[move][row][i];
        double lgCurProb = std::log(curProb);
        if (!std::isfinite(lgCurProb)) continue;
        if (moment == static_cast<uint8_t>(MomentType::FIRST))
            expectedLL += curProb * lgCurProb;
        else if (moment == static_cast<uint8_t>(MomentType::SECOND))
            expectedLL += curProb * (lgCurProb * lgCurProb);
    }
    return expectedLL;
}

PwSnrModel::PwSnrModel(const PwSnrModelCreator* params, const SNR& snr) : params_{params}, snr_(snr)
{
    for (size_t ctx = 0; ctx < CONTEXT_NUMBER; ++ctx) {
        const uint8_t bp = ctx & 3;  // (equivalent to % 4)
        const double snr1 = clip(snr_[bp], params_->snrRanges_[bp]), snr2 = snr1 * snr1,
                     snr3 = snr2 * snr1;
        double sum = 1.0;

        // cached transitions
        ctxTrans_[ctx][0] = 1.0;
        for (size_t j = 0; j < 3; ++j) {
            double xb = params_->transitionParams_[ctx][j][0] +
                        params_->transitionParams_[ctx][j][1] * snr1 +
                        params_->transitionParams_[ctx][j][2] * snr2 +
                        params_->transitionParams_[ctx][j][3] * snr3;
            xb = std::exp(xb);
            ctxTrans_[ctx][j + 1] = xb;
            sum += xb;
        }
        for (size_t j = 0; j < 4; ++j)
            ctxTrans_[ctx][j] /= sum;

        // cached expectations
        for (size_t move = 0; move < 3; ++move)
            for (size_t moment = 0; moment < 2; ++moment)
                cachedEmissionExpectations_[ctx][move][moment] =
                    CalculateExpectedLLForEmission(move, ctx, moment);
    }
}

std::unique_ptr<AbstractRecursor> PwSnrModel::CreateRecursor(const MappedRead& mr,
                                                             double scoreDiff) const
{
    const double counterWeight = CounterWeight(
        [this](size_t ctx, MoveType m) { return ctxTrans_[ctx][static_cast<uint8_t>(m)]; },
        [this](size_t ctx, MoveType m) {
            double r = 0.0;
            for (size_t o = 0; o < OUTCOME_NUMBER; ++o) {
                const double p = params_->emissionPmf_[static_cast<uint8_t>(m)][ctx][o];
                if (p > 0.0) r += p * std::log(p);
            }
            return r;
        },
        CONTEXT_NUMBER);

    return std::make_unique<PwSnrRecursor>(mr, scoreDiff, counterWeight, params_);
}

std::vector<TemplatePosition> PwSnrModel::Populate(const std::string& tpl) const
{
    auto rowFetcher = [this](const NCBI2na prev, const NCBI2na curr) -> const double(&)[4]
    {
        const auto row = EncodeContext16(prev, curr);
        const double(&params)[4] = ctxTrans_[row];
        return params;
    };
    return AbstractPopulater(tpl, rowFetcher);
}

double PwSnrModel::ExpectedLLForEmission(const MoveType move, const AlleleRep& prev,
                                         const AlleleRep& curr, const MomentType moment) const
{
    auto cachedEmissionVisitor = [this](const MoveType move, const NCBI2na prev, const NCBI2na curr,
                                        const MomentType moment) -> double {
        const auto row = EncodeContext16(prev, curr);
        return cachedEmissionExpectations_[row][static_cast<uint8_t>(move)]
                                          [static_cast<uint8_t>(moment)];
    };
    return AbstractExpectedLLForEmission(move, prev, curr, moment, cachedEmissionVisitor);
}

PwSnrRecursor::PwSnrRecursor(const MappedRead& mr, double scoreDiff, double counterWeight,
                             const PwSnrModelCreator* params)
    : Recursor(mr, scoreDiff)
    , params_{params}
    , counterWeight_{counterWeight}
    , nLgCounterWeight_{-std::log(counterWeight_)}
{
}

std::vector<uint8_t> PwSnrRecursor::EncodeRead(const MappedRead& read)
{
    std::vector<uint8_t> result;
    result.reserve(read.Length());

    for (size_t i = 0; i < read.Length(); ++i) {
        result.emplace_back(EncodeBase(read.Seq[i], read.PulseWidth[i]));
    }

    return result;
}

double PwSnrRecursor::EmissionPr(const MoveType move, const uint8_t emission, const AlleleRep& prev,
                                 const AlleleRep& curr) const
{
    return AbstractEmissionPr(params_->emissionPmf_, move, emission, prev, curr) * counterWeight_;
}

double PwSnrRecursor::UndoCounterWeights(const size_t nEmissions) const
{
    return nLgCounterWeight_ * nEmissions;
}

PwSnrModelCreator::PwSnrModelCreator(const boost::property_tree::ptree& pt)
{
    try {
        ReadMatrix<4, 2>(snrRanges_, pt.get_child("SnrRanges"));
        ReadMatrix<3, CONTEXT_NUMBER, OUTCOME_NUMBER>(emissionPmf_,
                                                      pt.get_child("EmissionParameters"));
        ReadMatrix<CONTEXT_NUMBER, 3, 4>(transitionParams_, pt.get_child("TransitionParameters"));
    } catch (std::invalid_argument& e) {
        throw MalformedModelFile();
    } catch (boost::property_tree::ptree_error&) {
        throw MalformedModelFile();
    }
}

class PwSnrInitializeModel
{
public:
    PwSnrInitializeModel(const PwSnrModel& model) : model_(model) {}

    inline std::pair<Data::SNR, std::vector<TemplatePosition>> operator()(
        std::default_random_engine* const rng, const std::string& tpl)
    {
        Data::SNR snrs{0, 0, 0, 0};
        for (uint8_t i = 0; i < 4; ++i) {
            snrs[i] = std::uniform_real_distribution<double>{
                model_.params_->snrRanges_[i][0], model_.params_->snrRanges_[i][1]}(*rng);
        }

        std::vector<TemplatePosition> transModel = model_.Populate(tpl);

        return {snrs, transModel};
    }

private:
    const PwSnrModel& model_;
};

class PwSnrGenerateReadData
{
public:
    PwSnrGenerateReadData(const PwSnrModelCreator& params) : params_(params) {}

    BaseData operator()(std::default_random_engine* const rng, const MoveType state,
                        const AlleleRep& prev, const AlleleRep& curr)
    {
        // distribution is arbitrary at the moment, as
        // IPD is not a covariate of the consensus HMM
        std::uniform_int_distribution<uint8_t> ipdDistrib(1, 5);

        std::array<double, OUTCOME_NUMBER> emissionDist;
        for (size_t i = 0; i < OUTCOME_NUMBER; ++i) {
            emissionDist[i] = AbstractEmissionPr(params_.emissionPmf_, state, i, prev, curr);
        }

        std::discrete_distribution<uint8_t> outcomeDistrib(emissionDist.cbegin(),
                                                           emissionDist.cend());

        const uint8_t event = outcomeDistrib(*rng);
        const std::pair<char, uint8_t> outcome = DecodeEmission(event);

        return {outcome.first, outcome.second, ipdDistrib(*rng)};
    }

private:
    const PwSnrModelCreator& params_;
};

std::pair<Data::Read, std::vector<MoveType>> PwSnrModel::SimulateRead(
    std::default_random_engine* const rng, const std::string& tpl,
    const std::string& readname) const
{
    const PwSnrInitializeModel init(*this);
    const PwSnrGenerateReadData generateData(*params_);

    return SimulateReadImpl(rng, tpl, readname, init, generateData);
}

}  // namespace anonymous
}  // namespace PwSnr

REGISTER_MODELFORM_IMPL(PwSnr)

}  // namespace Consensus
}  // namespace PacBio
