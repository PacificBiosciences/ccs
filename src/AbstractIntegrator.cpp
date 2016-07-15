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
#include <pacbio/consensus/Sequence.h>

#include "ModelFactory.h"

namespace PacBio {
namespace Consensus {

namespace {

constexpr double NEG_DBL_INF = -std::numeric_limits<double>::infinity();
constexpr int NEG_INT_INF = -std::numeric_limits<int>::infinity();
constexpr float NEG_FLOAT_INF = -std::numeric_limits<float>::infinity();

}  // anonymous namespace

std::set<std::string> SupportedChemistries() { return ModelFactory::SupportedChemistries(); }
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
State AbstractIntegrator::AddRead(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& read)
{
    // TODO(atoepfer) Why don't we add those reads and tag them as TEMPLATE_TOO_SMALL
    //                and effectively keep book about them? This logic should be
    //                in the Evaluator
    if (read.TemplateEnd <= read.TemplateStart) throw std::invalid_argument("template span < 2!");

    if (read.Length() < 2) throw std::invalid_argument("read span < 2!");

    evals_.emplace_back(Evaluator(std::move(tpl), read, cfg_.MinZScore, cfg_.ScoreDiff));
    return evals_.back().Status();
}

double AbstractIntegrator::LL(const Mutation& fwdMut) { return AccumulateNoInf(LLs(fwdMut)); }
double AbstractIntegrator::LL() const { return AccumulateNoInf(LLs()); }
std::vector<double> AbstractIntegrator::LLs(const Mutation& fwdMut)
{
    const Mutation revMut(ReverseComplement(fwdMut));

    const auto functor = [&fwdMut, &revMut](Evaluator& eval) {
        switch (eval.Strand()) {
            case StrandType::FORWARD:
                return eval.LL(fwdMut);
            case StrandType::REVERSE:
                return eval.LL(revMut);
            case StrandType::UNMAPPED:
                return NEG_DBL_INF;
            default:
                throw std::runtime_error("Unknown StrandType");
        }
    };

    return TransformEvaluators<double>(functor);
}

std::vector<double> AbstractIntegrator::LLs() const
{
    const auto functor = [](const Evaluator& eval) { return eval.LL(); };
    return TransformEvaluators<double>(functor);
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
std::vector<float> AbstractIntegrator::AlphaPopulated() const
{
    return TransformEvaluators<float>([](const Evaluator& eval) { return eval.AlphaPopulated(); });
}

float AbstractIntegrator::MaxAlphaPopulated() const { return MaxElement<float>(AlphaPopulated()); }
std::vector<float> AbstractIntegrator::BetaPopulated() const
{
    return TransformEvaluators<float>([](const Evaluator& eval) { return eval.BetaPopulated(); });
}

float AbstractIntegrator::MaxBetaPopulated() const { return MaxElement<float>(BetaPopulated()); }
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

std::vector<State> AbstractIntegrator::States() const
{
    return TransformEvaluators<State>([](const Evaluator& eval) { return eval.Status(); });
}

std::vector<StrandType> AbstractIntegrator::StrandTypes() const
{
    return TransformEvaluators<StrandType>([](const Evaluator& eval) { return eval.Strand(); });
}

Mutation AbstractIntegrator::ReverseComplement(const Mutation& mut) const
{
    return Mutation(mut.Type, TemplateLength() - mut.End(), Complement(mut.Base));
}

}  // namespace Consensus
}  // namespace PacBio
