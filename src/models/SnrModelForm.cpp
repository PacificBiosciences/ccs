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

template <typename T>
inline T clip(const T val, const T (&range)[2])
{
    return std::max(range[0], std::min(val, range[1]));
}

constexpr size_t CONTEXT_NUMBER = 8;

// fwd decl
class SnrModelCreator;

class SnrModel : public ModelConfig
{
public:
    SnrModel(const SnrModelCreator* params, const SNR& snr);
    std::unique_ptr<AbstractRecursor> CreateRecursor(const MappedRead& mr, double scoreDiff) const;
    std::vector<TemplatePosition> Populate(const std::string& tpl) const;
    double ExpectedLLForEmission(MoveType move, uint8_t prev, uint8_t curr,
                                 MomentType moment) const;

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
    double EmissionPr(MoveType move, uint8_t emission, uint8_t prev, uint8_t curr) const;
    double UndoCounterWeights(size_t nEmissions) const;

private:
    const SnrModelCreator* params_;
    double counterWeight_;
    double nLgCounterWeight_;
};

class SnrModelCreator : public ModelCreator
{
    REGISTER_MODELFORM(SnrModelCreator);
    friend class SnrModel;
    friend class SnrRecursor;

public:
    static ModelForm Form() { return ModelForm::SNR; }
    SnrModelCreator(const boost::property_tree::ptree& pt);
    virtual std::unique_ptr<ModelConfig> Create(const SNR& snr) const
    {
        return std::unique_ptr<ModelConfig>(new SnrModel(this, snr));
    };

private:
    double snrRanges_[4][2];
    double emissionPmf_[3][2];
    double transitionParams_[CONTEXT_NUMBER][3][4];
    double substitutionRate_;
};

REGISTER_MODELFORM_IMPL(SnrModelCreator);

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

    return std::unique_ptr<AbstractRecursor>(
        new SnrRecursor(mr, scoreDiff, counterWeight, params_));
}

std::vector<TemplatePosition> SnrModel::Populate(const std::string& tpl) const
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
        const uint8_t row = ((prev == curr) << 2) | curr;
        const auto& params = ctxTrans_[row];
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

double SnrModel::ExpectedLLForEmission(const MoveType move, const uint8_t prev, const uint8_t curr,
                                       const MomentType moment) const
{
    const double lgThird = -std::log(3.0);
    if (move == MoveType::MATCH) {
        const double probMismatch = params_->substitutionRate_;
        const double probMatch = 1.0 - probMismatch;
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
        const uint8_t em = detail::TranslationTable[static_cast<uint8_t>(bp)];
        if (em > 3) throw std::invalid_argument("invalid character in read!");
        result.emplace_back(em);
    }

    return result;
}

double SnrRecursor::EmissionPr(MoveType move, uint8_t emission, uint8_t prev, uint8_t curr) const
{
    assert(move != MoveType::DELETION);

    return params_->emissionPmf_[static_cast<uint8_t>(move)][curr != emission] * counterWeight_;
}

double SnrRecursor::UndoCounterWeights(const size_t nEmissions) const
{
    return nLgCounterWeight_ * nEmissions;
}

SnrModelCreator::SnrModelCreator(const boost::property_tree::ptree& pt)
    : emissionPmf_{{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0 / 3.0}}
{
    try {
        ReadMatrix<4, 2>(snrRanges_, pt.get_child("SnrRanges"));
        ReadMatrix<CONTEXT_NUMBER, 3, 4>(transitionParams_, pt.get_child("TransitionParameters"));
        substitutionRate_ = pt.get<double>("SubstitutionRate");
        emissionPmf_[0][0] = 1.0 - substitutionRate_;
        emissionPmf_[0][1] = substitutionRate_ / 3.0;
    } catch (std::invalid_argument& e) {
        throw MalformedModelFile();
    } catch (boost::property_tree::ptree_error&) {
        throw MalformedModelFile();
    }
}

}  // namespace anonymous
}  // namespace Consensus
}  // namespace PacBio
