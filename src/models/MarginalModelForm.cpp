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
namespace Marginal {
namespace {

using MalformedModelFile = PacBio::Exception::MalformedModelFile;

static constexpr const size_t CONTEXT_NUMBER = 8;
static constexpr const size_t OUTCOME_NUMBER = 4;

static constexpr const double snrRanges[2][4] = {
    {0, 0, 0, 0},     // minimum
    {20, 19, 20, 20}  // maximum
};

// fwd decl
class MarginalModelCreator;

class MarginalModel : public ModelConfig
{
public:
    MarginalModel(const MarginalModelCreator* params, const SNR& snr);
    std::unique_ptr<AbstractRecursor> CreateRecursor(const MappedRead& mr,
                                                     double scoreDiff) const override;
    std::vector<TemplatePosition> Populate(const std::string& tpl) const override;
    std::pair<Data::Read, std::vector<MoveType>> SimulateRead(
        std::default_random_engine* const rng, const std::string& tpl,
        const std::string& readname) const override;
    double ExpectedLLForEmission(MoveType move, const AlleleRep& prev, const AlleleRep& curr,
                                 MomentType moment) const override;

private:
    const MarginalModelCreator* params_;
};

class MarginalRecursor : public Recursor<MarginalRecursor>
{
public:
    MarginalRecursor(const MappedRead& mr, double scoreDiff, double counterWeight,
                     const MarginalModelCreator* params);

    static std::vector<uint8_t> EncodeRead(const MappedRead& read);
    double EmissionPr(MoveType move, uint8_t emission, const AlleleRep& prev,
                      const AlleleRep& curr) const;
    double UndoCounterWeights(size_t nEmissions) const override;

private:
    const MarginalModelCreator* params_;
    double counterWeight_;
    double nLgCounterWeight_;
};

class MarginalModelCreator : public ModelCreator
{
    friend class MarginalModel;
    friend class MarginalRecursor;
    friend class MarginalModelGenerateReadData;

public:
    static ModelForm Form() { return ModelForm::MARGINAL; }
    MarginalModelCreator(const boost::property_tree::ptree& pt);
    std::unique_ptr<ModelConfig> Create(const SNR& snr) const override
    {
        return std::make_unique<MarginalModel>(this, snr);
    };

private:
    double emissionPmf_[3][CONTEXT_NUMBER][OUTCOME_NUMBER];
    double transitionPmf_[CONTEXT_NUMBER][4];
};

MarginalModel::MarginalModel(const MarginalModelCreator* params, const SNR& snr) : params_{params}
{
}

std::unique_ptr<AbstractRecursor> MarginalModel::CreateRecursor(const MappedRead& mr,
                                                                double scoreDiff) const
{
    const double counterWeight = CounterWeight(
        [this](size_t ctx, MoveType m) {
            return params_->transitionPmf_[ctx][static_cast<uint8_t>(m)];
        },
        [this](size_t ctx, MoveType m) {
            double r = 0.0;
            for (size_t o = 0; o < OUTCOME_NUMBER; ++o) {
                const double p = params_->emissionPmf_[static_cast<uint8_t>(m)][ctx][o];
                if (p > 0.0) r += p * std::log(p);
            }
            return r;
        },
        CONTEXT_NUMBER);

    return std::make_unique<MarginalRecursor>(mr, scoreDiff, counterWeight, params_);
}

std::vector<TemplatePosition> MarginalModel::Populate(const std::string& tpl) const
{
    auto rowFetcher = [this](const NCBI2na prev, const NCBI2na curr) -> const double(&)[4]
    {
        const auto ctx = EncodeContext8(prev, curr);
        const double(&params)[4] = params_->transitionPmf_[ctx];
        return params;
    };
    return AbstractPopulater(tpl, rowFetcher);
}

double MarginalModel::ExpectedLLForEmission(const MoveType move, const AlleleRep& prev,
                                            const AlleleRep& curr, const MomentType moment) const
{
    auto cachedEmissionVisitor = [this](const MoveType move, const NCBI2na prev, const NCBI2na curr,
                                        const MomentType moment) -> double {
        const auto ctx = EncodeContext8(prev, curr);
        double expectedLL = 0;
        for (size_t i = 0; i < OUTCOME_NUMBER; i++) {
            double curProb = params_->emissionPmf_[static_cast<uint8_t>(move)][ctx][i];
            double lgCurProb = std::log(curProb);
            if (!std::isfinite(lgCurProb)) continue;
            if (moment == MomentType::FIRST)
                expectedLL += curProb * lgCurProb;
            else if (moment == MomentType::SECOND)
                expectedLL += curProb * (lgCurProb * lgCurProb);
        }
        return expectedLL;
    };
    return AbstractExpectedLLForEmission(move, prev, curr, moment, cachedEmissionVisitor);
}

MarginalRecursor::MarginalRecursor(const MappedRead& mr, double scoreDiff, double counterWeight,
                                   const MarginalModelCreator* params)
    : Recursor(mr, scoreDiff)
    , params_{params}
    , counterWeight_{counterWeight}
    , nLgCounterWeight_{-std::log(counterWeight_)}
{
}

std::vector<uint8_t> MarginalRecursor::EncodeRead(const MappedRead& read)
{
    std::vector<uint8_t> result;
    result.reserve(read.Length());

    for (const char bp : read.Seq) {
        result.emplace_back(EncodeBase(bp));
    }

    return result;
}

double MarginalRecursor::EmissionPr(const MoveType move, const uint8_t emission,
                                    const AlleleRep& prev, const AlleleRep& curr) const
{
    return AbstractEmissionPr(params_->emissionPmf_, move, emission, prev, curr) * counterWeight_;
}

double MarginalRecursor::UndoCounterWeights(const size_t nEmissions) const
{
    return nLgCounterWeight_ * nEmissions;
}

MarginalModelCreator::MarginalModelCreator(const boost::property_tree::ptree& pt)
{
    try {
        ReadMatrix<3, CONTEXT_NUMBER, OUTCOME_NUMBER>(emissionPmf_,
                                                      pt.get_child("EmissionParameters"));
        ReadMatrix<CONTEXT_NUMBER, 4>(transitionPmf_, pt.get_child("TransitionParameters"));
    } catch (std::invalid_argument& e) {
        throw MalformedModelFile();
    } catch (boost::property_tree::ptree_error&) {
        throw MalformedModelFile();
    }
}

class MarginalModelInitializeModel
{
public:
    MarginalModelInitializeModel(const MarginalModel& model) : model_(model) {}

    inline std::pair<Data::SNR, std::vector<TemplatePosition>> operator()(
        std::default_random_engine* const rng, const std::string& tpl)
    {
        Data::SNR snrs{0, 0, 0, 0};
        for (uint8_t i = 0; i < 4; ++i) {
            snrs[i] =
                std::uniform_real_distribution<double>{snrRanges[0][i], snrRanges[1][i]}(*rng);
        }

        std::vector<TemplatePosition> transModel = model_.Populate(tpl);

        return {snrs, transModel};
    }

private:
    const MarginalModel& model_;
};

class MarginalModelGenerateReadData
{
public:
    MarginalModelGenerateReadData(const MarginalModelCreator& params) : params_(params) {}

    BaseData operator()(std::default_random_engine* const rng, const MoveType state,
                        const AlleleRep& prev, const AlleleRep& curr)
    {
        // distribution is arbitrary at the moment, as
        // PW and IPD are not a covariates of the consensus HMM
        std::uniform_int_distribution<uint8_t> pwDistrib{1, 3};
        std::uniform_int_distribution<uint8_t> ipdDistrib{1, 5};

        std::array<double, 4> baseDist;
        for (size_t i = 0; i < 4; ++i) {
            baseDist[i] = AbstractEmissionPr(params_.emissionPmf_, state, i, prev, curr);
        }

        std::discrete_distribution<uint8_t> baseDistrib(baseDist.cbegin(), baseDist.cend());

        const char newBase = Data::detail::NCBI2naToASCIIImpl(baseDistrib(*rng));
        const uint8_t newPw = pwDistrib(*rng);
        const uint8_t newIpd = ipdDistrib(*rng);

        return {newBase, newPw, newIpd};
    }

private:
    const MarginalModelCreator& params_;
};

std::pair<Data::Read, std::vector<MoveType>> MarginalModel::SimulateRead(
    std::default_random_engine* const rng, const std::string& tpl,
    const std::string& readname) const
{
    const MarginalModelInitializeModel init(*this);
    const MarginalModelGenerateReadData generateData(*params_);

    return SimulateReadImpl(rng, tpl, readname, init, generateData);
}

}  // namespace anonymous
}  // namespace Marginal

REGISTER_MODELFORM_IMPL(Marginal)

}  // namespace Consensus
}  // namespace PacBio
