// Copyright (c) 2011-2016, Pacific Biosciences of California, Inc.
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

#include <cassert>
#include <cmath>
#include <exception>
#include <limits>
#include <numeric>
#include <utility>

#include <pacbio/consensus/AbstractIntegrator.h>
#include <pacbio/consensus/AbstractMatrix.h>
#include <pacbio/exception/InvalidEvaluatorException.h>

#include <pacbio/data/Sequence.h>

#include "ModelFactory.h"

namespace PacBio {
namespace Consensus {

using namespace PacBio::Data;
using namespace PacBio::Exception;

namespace {

constexpr double NEG_DBL_INF = -std::numeric_limits<double>::infinity();
constexpr int NEG_INT_INF = -std::numeric_limits<int>::infinity();
constexpr float NEG_FLOAT_INF = -std::numeric_limits<float>::infinity();

}  // anonymous namespace

IntegratorConfig::IntegratorConfig(const double minZScore, const double scoreDiff)
    : MinZScore{minZScore}, ScoreDiff{scoreDiff}
{
    if (ScoreDiff < 0) {
        throw std::runtime_error("Score diff must be > 0");
    }
}

AbstractIntegrator::AbstractIntegrator(const IntegratorConfig& cfg) : cfg_{cfg} {}

AbstractIntegrator::AbstractIntegrator(AbstractIntegrator&& ai)
    : cfg_{ai.cfg_}, evals_{std::move(ai.evals_)}
{
}

AbstractIntegrator::~AbstractIntegrator() {}

Data::State AbstractIntegrator::AddRead(std::unique_ptr<AbstractTemplate>&& tpl,
                                        const Data::MappedRead& read)
{
    // TODO(atoepfer) Why don't we add those reads and tag them as TEMPLATE_TOO_SMALL
    //                and effectively keep book about them? This logic should be
    //                in the Evaluator
    if (read.TemplateEnd <= read.TemplateStart) throw std::invalid_argument("template span < 2!");

    if (read.Length() < 2) throw std::invalid_argument("read span < 2!");

    evals_.emplace_back(Evaluator(std::move(tpl), read, cfg_.MinZScore, cfg_.ScoreDiff));
    return evals_.back().Status();
}

double AbstractIntegrator::LL(const Mutation& fwdMut)
{
    const auto lls = LLs(fwdMut);
    return std::accumulate(lls.cbegin(), lls.cend(), 0.0);
}

double AbstractIntegrator::LL() const
{
    const auto functor = [](const Evaluator& eval) { return eval.IsValid() ? eval.LL() : 0; };
    const auto lls = TransformEvaluators<double>(functor);
    return std::accumulate(lls.cbegin(), lls.cend(), 0.0);
}

std::vector<double> AbstractIntegrator::LLs(const Mutation& fwdMut)
{
    const Mutation revMut(ReverseComplement(fwdMut));

    // Compute individual LLs of each Evaluator
    std::vector<double> lls;
    lls.reserve(evals_.size());
    for (auto& e : evals_) {
        // Ignore invalid Evaluators
        if (!e.IsValid()) continue;

        double ll;

        switch (e.Strand()) {
            case StrandType::FORWARD:
                ll = e.LL(fwdMut);
                break;
            case StrandType::REVERSE:
                ll = e.LL(revMut);
                break;
            case StrandType::UNMAPPED:
                // unmapped Evaluators should not be used
                throw InvalidEvaluatorException("Unmapped read in mutation testing");
            default:
                throw std::runtime_error("Unknown StrandType");
        }

        lls.emplace_back(ll);
    }

    return lls;
}

std::vector<std::string> AbstractIntegrator::ReadNames() const
{
    return TransformEvaluators<std::string>([](const Evaluator& eval) { return eval.ReadName(); });
}

std::vector<int> AbstractIntegrator::NumFlipFlops() const
{
    return TransformEvaluators<int>([](const Evaluator& eval) { return eval.NumFlipFlops(); });
}

int AbstractIntegrator::MaxNumFlipFlops() const { return MaxElement<int>(NumFlipFlops()); }

float AbstractIntegrator::MaxAlphaPopulated() const
{
    const auto functor = [](const Evaluator& eval) {
        return (eval.IsValid() ? eval.Alpha().UsedEntriesRatio() : NEG_FLOAT_INF);
    };
    auto alphaPopulated = TransformEvaluators<float>(functor);
    return MaxElement<float>(alphaPopulated);
}

float AbstractIntegrator::MaxBetaPopulated() const
{
    const auto functor = [](const Evaluator& eval) {
        return (eval.IsValid() ? eval.Beta().UsedEntriesRatio() : NEG_FLOAT_INF);
    };
    auto betaPopulated = TransformEvaluators<float>(functor);
    return MaxElement<float>(betaPopulated);
}

double AbstractIntegrator::AvgZScore() const
{
    double mean = 0.0, var = 0.0;
    size_t n = 0;
    for (const auto& eval : evals_) {
        if (eval) {
            double m, v;
            std::tie(m, v) = eval.NormalParameters();
            mean += m;
            var += v;
            ++n;
        }
    }
    return (LL() / n - mean / n) / std::sqrt(var / n);
}

std::vector<double> AbstractIntegrator::ZScores() const
{
    return TransformEvaluators<double>([](const Evaluator& eval) { return eval.ZScore(); });
}

std::vector<std::pair<double, double>> AbstractIntegrator::NormalParameters() const
{
    return TransformEvaluators<std::pair<double, double>>(
        [](const Evaluator& eval) { return eval.NormalParameters(); });
}

std::vector<Data::State> AbstractIntegrator::States() const
{
    return TransformEvaluators<Data::State>([](const Evaluator& eval) { return eval.Status(); });
}

std::vector<StrandType> AbstractIntegrator::StrandTypes() const
{
    return TransformEvaluators<StrandType>([](const Evaluator& eval) { return eval.Strand(); });
}

const Evaluator& AbstractIntegrator::GetEvaluator(size_t idx) const { return evals_[idx]; }

const AbstractMatrix& AbstractIntegrator::Alpha(size_t idx) const { return evals_[idx].Alpha(); }

const AbstractMatrix& AbstractIntegrator::Beta(size_t idx) const { return evals_[idx].Beta(); }

Mutation AbstractIntegrator::ReverseComplement(const Mutation& mut) const
{
    return Mutation(mut.Type, TemplateLength() - mut.End(), Complement(mut.Base));
}

}  // namespace Consensus
}  // namespace PacBio
