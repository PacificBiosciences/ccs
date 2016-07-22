// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
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

// Author: Lance Hepler

#include <cassert>
#include <cmath>
#include <memory>
#include <stdexcept>

#include <pacbio/consensus/ModelConfig.h>
#include <pacbio/data/Read.h>
#include <pacbio/exception/ModelError.h>

#include "../JsonHelpers.h"
#include "../ModelFactory.h"
#include "../ModelFormFactory.h"
#include "../Recursor.h"

using namespace PacBio::Data;

namespace PacBio {
namespace Consensus {
namespace {

using MalformedModelFile = PacBio::Exception::MalformedModelFile;

template <typename T>
inline T clip(const T val, const T range[2])
{
    return std::max(range[0], std::min(val, range[1]));
}

constexpr size_t OUTCOME_NUMBER = 12;
constexpr size_t CONTEXT_NUMBER = 16;

// fwd decl
class PwSnrAModelCreator;

class PwSnrAModel : public ModelConfig
{
public:
    PwSnrAModel(const PwSnrAModelCreator* params, const SNR& snr);
    std::unique_ptr<AbstractRecursor> CreateRecursor(std::unique_ptr<AbstractTemplate>&& tpl,
                                                     const MappedRead& mr, double scoreDiff) const;
    std::vector<TemplatePosition> Populate(const std::string& tpl) const;
    double ExpectedLLForEmission(MoveType move, uint8_t prev, uint8_t curr,
                                 MomentType moment) const;

private:
    double CalculateExpectedLLForEmission(const size_t move, const uint8_t row,
                                          const size_t moment) const;

    const PwSnrAModelCreator* params_;
    SNR snr_;
    double ctxTrans_[CONTEXT_NUMBER][4];
    double cachedEmissionExpectations_[CONTEXT_NUMBER][3][2];
    double counterWeight_;
};

class PwSnrARecursor : public Recursor<PwSnrARecursor>
{
public:
    PwSnrARecursor(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr, double scoreDiff,
                   const PwSnrAModelCreator* params, double counterWeight);

    static std::vector<uint8_t> EncodeRead(const MappedRead& read);
    double EmissionPr(MoveType move, uint8_t emission, uint8_t prev, uint8_t curr) const;
    double UndoCounterWeights(size_t nEmissions) const;

private:
    const PwSnrAModelCreator* params_;
    double counterWeight_;
    double nLgCounterWeight_;
};

class PwSnrAModelCreator : public ModelCreator
{
    REGISTER_MODELFORM(PwSnrAModelCreator);
    friend class PwSnrAModel;
    friend class PwSnrARecursor;

public:
    static std::string Name() { return "PwSnrA"; }
    PwSnrAModelCreator(const boost::property_tree::ptree& pt);
    virtual std::unique_ptr<ModelConfig> Create(const SNR& snr) const
    {
        return std::unique_ptr<ModelConfig>(new PwSnrAModel(this, snr));
    };

private:
    double snrRanges_[2];
    double emissionPmf_[3][CONTEXT_NUMBER][OUTCOME_NUMBER];
    double transitionParams_[CONTEXT_NUMBER][3][4];
    double counterWeight_;
};

REGISTER_MODELFORM_IMPL(PwSnrAModelCreator);

inline double PwSnrAModel::CalculateExpectedLLForEmission(const size_t move, const uint8_t row,
                                                          const size_t moment) const
{
    double expectedLL = 0;
    for (size_t i = 0; i < OUTCOME_NUMBER; i++) {
        double curProb = params_->emissionPmf_[move][row][i];
        double lgCurProb = std::log(curProb);
        if (moment == static_cast<uint8_t>(MomentType::FIRST))
            expectedLL += curProb * lgCurProb;
        else if (moment == static_cast<uint8_t>(MomentType::SECOND))
            expectedLL += curProb * (lgCurProb * lgCurProb);
    }
    return expectedLL;
}

PwSnrAModel::PwSnrAModel(const PwSnrAModelCreator* params, const SNR& snr)
    : params_{params}, snr_(snr)
{
    const double snr1 = clip(snr_.A, params_->snrRanges_), snr2 = snr1 * snr1, snr3 = snr2 * snr1;
    for (size_t ctx = 0; ctx < CONTEXT_NUMBER; ++ctx) {
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

    counterWeight_ = params_->counterWeight_;
}

std::unique_ptr<AbstractRecursor> PwSnrAModel::CreateRecursor(
    std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr, double scoreDiff) const
{
    return std::unique_ptr<AbstractRecursor>(
        new PwSnrARecursor(std::forward<std::unique_ptr<AbstractTemplate>>(tpl), mr, scoreDiff,
                           params_, counterWeight_));
}

std::vector<TemplatePosition> PwSnrAModel::Populate(const std::string& tpl) const
{
    std::vector<TemplatePosition> result;

    if (tpl.empty()) return result;

    result.reserve(tpl.size());

    // calculate transition probabilities
    uint8_t prev = detail::TranslationTable[static_cast<uint8_t>(tpl[0])];
    if (prev > 3) throw std::invalid_argument("invalid character in template!");

    for (size_t i = 1; i < tpl.size(); ++i) {
        const uint8_t curr = detail::TranslationTable[static_cast<uint8_t>(tpl[i])];
        if (curr > 3) throw std::invalid_argument("invalid character in template!");
        const auto row = (prev << 2) | curr;
        const auto params = ctxTrans_[row];
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

double PwSnrAModel::ExpectedLLForEmission(const MoveType move, const uint8_t prev,
                                          const uint8_t curr, const MomentType moment) const
{
    const size_t row = (prev << 2) | curr;
    return cachedEmissionExpectations_[row][static_cast<uint8_t>(move)]
                                      [static_cast<uint8_t>(moment)];
}

PwSnrARecursor::PwSnrARecursor(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr,
                               double scoreDiff, const PwSnrAModelCreator* params,
                               double counterWeight)
    : Recursor(std::forward<std::unique_ptr<AbstractTemplate>>(tpl), mr, scoreDiff)
    , params_{params}
    , counterWeight_{counterWeight}
    , nLgCounterWeight_{-std::log(counterWeight_)}
{
}

std::vector<uint8_t> PwSnrARecursor::EncodeRead(const MappedRead& read)
{
    std::vector<uint8_t> result;

    result.reserve(read.Length());

    for (size_t i = 0; i < read.Length(); ++i) {
        if (read.PulseWidth[i] < 1U) throw std::runtime_error("invalid PulseWidth in read!");
        const uint8_t pw = std::min(2, read.PulseWidth[i] - 1);
        const uint8_t bp = detail::TranslationTable[static_cast<uint8_t>(read.Seq[i])];
        if (bp > 3) throw std::invalid_argument("invalid character in read!");
        const uint8_t em = (pw << 2) | bp;
        if (em > 11) throw std::runtime_error("read encoding error!");
        result.emplace_back(em);
    }

    return result;
}

double PwSnrARecursor::EmissionPr(MoveType move, uint8_t emission, uint8_t prev, uint8_t curr) const
{
    assert(move != MoveType::DELETION);
    const auto row = (prev << 2) | curr;
    return params_->emissionPmf_[static_cast<uint8_t>(move)][row][emission] * counterWeight_;
}

double PwSnrARecursor::UndoCounterWeights(const size_t nEmissions) const
{
    return nLgCounterWeight_ * nEmissions;
}

PwSnrAModelCreator::PwSnrAModelCreator(const boost::property_tree::ptree& pt)
{
    try {
        ReadMatrix<2>(snrRanges_, pt.get_child("SnrRanges"));
        ReadMatrix<3, CONTEXT_NUMBER, OUTCOME_NUMBER>(emissionPmf_,
                                                      pt.get_child("EmissionParameters"));
        ReadMatrix<CONTEXT_NUMBER, 3, 4>(transitionParams_, pt.get_child("TransitionParameters"));
        counterWeight_ = pt.get<double>("CounterWeight");
    } catch (std::invalid_argument& e) {
        throw MalformedModelFile();
    } catch (boost::property_tree::ptree_error) {
        throw MalformedModelFile();
    }
}

}  // namespace anonymous
}  // namespace Consensus
}  // namespace PacBio
