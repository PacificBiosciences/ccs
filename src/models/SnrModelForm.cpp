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
namespace Snr {
namespace {

using MalformedModelFile = PacBio::Exception::MalformedModelFile;

static constexpr const size_t CONTEXT_NUMBER = 8;

// fwd decl
class SnrModelCreator;

class SnrModel : public ModelConfig
{
public:
    SnrModel(const SnrModelCreator* params, const SNR& snr);
    std::unique_ptr<AbstractRecursor> CreateRecursor(const MappedRead& mr,
                                                     double scoreDiff) const override;
    std::vector<TemplatePosition> Populate(const std::string& tpl) const override;
    std::pair<Data::Read, std::vector<MoveType>> SimulateRead(
        std::default_random_engine* const rng, const std::string& tpl,
        const std::string& readname) const override;
    double ExpectedLLForEmission(MoveType move, const AlleleRep& prev, const AlleleRep& curr,
                                 MomentType moment) const override;
    friend class SnrInitializeModel;

private:
    const SnrModelCreator* params_;
    SNR snr_;
    double ctxTrans_[CONTEXT_NUMBER][4];
};

class SnrRecursor : public Recursor<SnrRecursor>
{
public:
    SnrRecursor(const MappedRead& mr, double scoreDiff, double counterWeight,
                const SnrModelCreator* params);

    static std::vector<uint8_t> EncodeRead(const MappedRead& read);
    double EmissionPr(MoveType move, uint8_t emission, const AlleleRep& prev,
                      const AlleleRep& curr) const;
    double UndoCounterWeights(size_t nEmissions) const override;

private:
    const SnrModelCreator* params_;
    double counterWeight_;
    double nLgCounterWeight_;
};

class SnrModelCreator : public ModelCreator
{
    friend class SnrModel;
    friend class SnrRecursor;
    friend class SnrInitializeModel;
    friend class SnrGenerateReadData;

public:
    static ModelForm Form() { return ModelForm::SNR; }
    SnrModelCreator(const boost::property_tree::ptree& pt);
    std::unique_ptr<ModelConfig> Create(const SNR& snr) const override
    {
        return std::make_unique<SnrModel>(this, snr);
    };

private:
    double snrRanges_[4][2];
    double emissionPmf_[3][1][2];
    double transitionParams_[CONTEXT_NUMBER][3][4];
    double substitutionRate_;
};

SnrModel::SnrModel(const SnrModelCreator* params, const SNR& snr) : params_{params}, snr_(snr)
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
    }
}

std::unique_ptr<AbstractRecursor> SnrModel::CreateRecursor(const MappedRead& mr,
                                                           double scoreDiff) const
{
    const double counterWeight = CounterWeight(
        [this](size_t ctx, MoveType m) { return ctxTrans_[ctx][static_cast<uint8_t>(m)]; },
        [this](size_t ctx, MoveType m) {
            const double kEps = params_->substitutionRate_;
            const double kInvEps = 1.0 - kEps;
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
        CONTEXT_NUMBER);

    return std::make_unique<SnrRecursor>(mr, scoreDiff, counterWeight, params_);
}

std::vector<TemplatePosition> SnrModel::Populate(const std::string& tpl) const
{
    auto rowFetcher = [this](const NCBI2na prev, const NCBI2na curr) -> const double(&)[4]
    {
        const auto row = EncodeContext8(prev, curr);
        const double(&params)[4] = ctxTrans_[row];
        return params;
    };
    return AbstractPopulater(tpl, rowFetcher);
}

double SnrModel::ExpectedLLForEmission(const MoveType move, const AlleleRep& prev,
                                       const AlleleRep& curr, const MomentType moment) const
{
    auto cachedEmissionVisitor = [this](const MoveType move, const NCBI2na prev, const NCBI2na curr,
                                        const MomentType moment) -> double {
        const double lgThird = -std::log(3.0);
        if (move == MoveType::MATCH) {
            const double probMismatch = params_->substitutionRate_;
            const double probMatch = 1.0 - probMismatch;
            const double lgMatch = std::log(probMatch);
            const double lgMismatch = lgThird + std::log(probMismatch);
            if (!std::isfinite(lgMatch) || !std::isfinite(lgMismatch)) return 0.0;
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

SnrRecursor::SnrRecursor(const MappedRead& mr, double scoreDiff, double counterWeight,
                         const SnrModelCreator* params)
    : Recursor(mr, scoreDiff)
    , params_{params}
    , counterWeight_{counterWeight}
    , nLgCounterWeight_{-std::log(counterWeight_)}
{
}

std::vector<uint8_t> SnrRecursor::EncodeRead(const MappedRead& read)
{
    std::vector<uint8_t> result;
    result.reserve(read.Length());

    for (const char bp : read.Seq) {
        result.emplace_back(EncodeBase(bp));
    }

    return result;
}

double SnrRecursor::EmissionPr(const MoveType move, const uint8_t emission, const AlleleRep& prev,
                               const AlleleRep& curr) const
{
    return AbstractEmissionPr(params_->emissionPmf_, move, emission, prev, curr) * counterWeight_;
}

double SnrRecursor::UndoCounterWeights(const size_t nEmissions) const
{
    return nLgCounterWeight_ * nEmissions;
}

SnrModelCreator::SnrModelCreator(const boost::property_tree::ptree& pt)
    : emissionPmf_{{{0.0, 0.0}}, {{1.0, 0.0}}, {{0.0, 1.0 / 3.0}}}
{
    try {
        ReadMatrix<4, 2>(snrRanges_, pt.get_child("SnrRanges"));
        ReadMatrix<CONTEXT_NUMBER, 3, 4>(transitionParams_, pt.get_child("TransitionParameters"));
        substitutionRate_ = pt.get<double>("SubstitutionRate");
        emissionPmf_[0][0][0] = 1.0 - substitutionRate_;
        emissionPmf_[0][0][1] = substitutionRate_ / 3.0;
    } catch (std::invalid_argument& e) {
        throw MalformedModelFile();
    } catch (boost::property_tree::ptree_error&) {
        throw MalformedModelFile();
    }
}

class SnrInitializeModel
{
public:
    SnrInitializeModel(const SnrModel& model) : model_(model) {}

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
    const SnrModel& model_;
};

class SnrGenerateReadData
{
public:
    SnrGenerateReadData(const SnrModelCreator& params) : params_(params) {}

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
    const SnrModelCreator& params_;
};

std::pair<Data::Read, std::vector<MoveType>> SnrModel::SimulateRead(
    std::default_random_engine* const rng, const std::string& tpl,
    const std::string& readname) const
{
    const SnrInitializeModel init(*this);
    const SnrGenerateReadData generateData(*params_);

    return SimulateReadImpl(rng, tpl, readname, init, generateData);
}

}  // namespace anonymous
}  // namespace Snr

REGISTER_MODELFORM_IMPL(Snr)

}  // namespace Consensus
}  // namespace PacBio
