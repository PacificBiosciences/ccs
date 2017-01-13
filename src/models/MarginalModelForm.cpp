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
#include "CounterWeight.h"

using namespace PacBio::Data;

namespace PacBio {
namespace Consensus {
namespace {

using MalformedModelFile = PacBio::Exception::MalformedModelFile;

constexpr size_t CONTEXT_NUMBER = 8;
constexpr size_t OUTCOME_NUMBER = 4;

// fwd decl
class MarginalModelCreator;

class MarginalModel : public ModelConfig
{
public:
    MarginalModel(const MarginalModelCreator* params, const SNR& snr);
    std::unique_ptr<AbstractRecursor> CreateRecursor(std::unique_ptr<AbstractTemplate>&& tpl,
                                                     const MappedRead& mr, double scoreDiff) const;
    std::vector<TemplatePosition> Populate(const std::string& tpl) const;
    double ExpectedLLForEmission(MoveType move, uint8_t prev, uint8_t curr,
                                 MomentType moment) const;

private:
    const MarginalModelCreator* params_;
};

class MarginalRecursor : public Recursor<MarginalRecursor>
{
public:
    MarginalRecursor(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr,
                     double scoreDiff, double counterWeight, const MarginalModelCreator* params);

    static std::vector<uint8_t> EncodeRead(const MappedRead& read);
    double EmissionPr(MoveType move, uint8_t emission, uint8_t prev, uint8_t curr) const;
    double UndoCounterWeights(size_t nEmissions) const;

private:
    const MarginalModelCreator* params_;
    double counterWeight_;
    double nLgCounterWeight_;
};

class MarginalModelCreator : public ModelCreator
{
    REGISTER_MODELFORM(MarginalModelCreator);
    friend class MarginalModel;
    friend class MarginalRecursor;

public:
    static std::string Name() { return "Marginal"; }
    MarginalModelCreator(const boost::property_tree::ptree& pt);
    virtual std::unique_ptr<ModelConfig> Create(const SNR& snr) const
    {
        return std::unique_ptr<ModelConfig>(new MarginalModel(this, snr));
    };

private:
    double emissionPmf_[3][CONTEXT_NUMBER][OUTCOME_NUMBER];
    double transitionPmf_[CONTEXT_NUMBER][4];
};

REGISTER_MODELFORM_IMPL(MarginalModelCreator);

MarginalModel::MarginalModel(const MarginalModelCreator* params, const SNR& snr) : params_{params}
{
}

std::unique_ptr<AbstractRecursor> MarginalModel::CreateRecursor(
    std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr, double scoreDiff) const
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

    return std::unique_ptr<AbstractRecursor>(
        new MarginalRecursor(std::forward<std::unique_ptr<AbstractTemplate>>(tpl), mr, scoreDiff,
                             counterWeight, params_));
}

std::vector<TemplatePosition> MarginalModel::Populate(const std::string& tpl) const
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
        const uint8_t ctx = ((prev == curr) << 2) | curr;
        const auto& params = params_->transitionPmf_[ctx];
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

double MarginalModel::ExpectedLLForEmission(const MoveType move, const uint8_t prev,
                                            const uint8_t curr, const MomentType moment) const
{
    const uint8_t ctx = ((prev == curr) << 2) | curr;
    double expectedLL = 0;
    for (size_t i = 0; i < OUTCOME_NUMBER; i++) {
        double curProb = params_->emissionPmf_[static_cast<uint8_t>(move)][ctx][i];
        double lgCurProb = std::log(curProb);
        if (moment == MomentType::FIRST)
            expectedLL += curProb * lgCurProb;
        else if (moment == MomentType::SECOND)
            expectedLL += curProb * (lgCurProb * lgCurProb);
    }
    return expectedLL;
}

MarginalRecursor::MarginalRecursor(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr,
                                   double scoreDiff, double counterWeight,
                                   const MarginalModelCreator* params)
    : Recursor(std::forward<std::unique_ptr<AbstractTemplate>>(tpl), mr, scoreDiff)
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
        const uint8_t em = detail::TranslationTable[static_cast<uint8_t>(bp)];
        if (em > 3) throw std::invalid_argument("invalid character in read!");
        result.emplace_back(em);
    }

    return result;
}

double MarginalRecursor::EmissionPr(MoveType move, uint8_t emission, uint8_t prev,
                                    uint8_t curr) const
{
    assert(move != MoveType::DELETION);
    const uint8_t ctx = ((prev == curr) << 2) | curr;
    return params_->emissionPmf_[static_cast<uint8_t>(move)][ctx][emission] * counterWeight_;
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
    } catch (boost::property_tree::ptree_error) {
        throw MalformedModelFile();
    }
}

}  // namespace anonymous
}  // namespace Consensus
}  // namespace PacBio
