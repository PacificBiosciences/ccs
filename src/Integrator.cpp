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

#include <pacbio/consensus/AbstractMatrix.h>
#include <pacbio/consensus/Integrator.h>
#include <pacbio/data/Sequence.h>
#include <pacbio/exception/InvalidEvaluatorException.h>

#include "Constants.h"
#include "ModelFactory.h"

using namespace PacBio::Data;
using namespace PacBio::Exception;

namespace PacBio {
namespace Consensus {

IntegratorConfig::IntegratorConfig(const double minZScore, const double scoreDiff)
    : MinZScore{minZScore}, ScoreDiff{scoreDiff}
{
    if (ScoreDiff < 0) {
        throw std::runtime_error("Score diff must be > 0");
    }
}

Integrator::Integrator(const std::string& tpl, const IntegratorConfig& cfg)
    : cfg_(cfg), fwdTpl_{tpl}, revTpl_{::PacBio::Data::ReverseComplement(tpl)}
{
}

Data::State Integrator::AddRead(std::unique_ptr<AbstractTemplate>&& tpl,
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

double Integrator::LL(const Mutation& fwdMut)
{
    const auto lls = LLs(fwdMut);
    return std::accumulate(lls.cbegin(), lls.cend(), 0.0);
}

double Integrator::LL() const
{
    const auto functor = [](const Evaluator& eval) { return eval.IsValid() ? eval.LL() : 0; };
    const auto lls = TransformEvaluators<double>(functor);
    return std::accumulate(lls.cbegin(), lls.cend(), 0.0);
}

std::vector<double> Integrator::LLs(const Mutation& fwdMut)
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

std::vector<double> Integrator::LLs() const
{
    const auto functor = [](const Evaluator& eval) { return eval.LL(); };
    return TransformEvaluators<double>(functor);
}

std::vector<std::string> Integrator::ReadNames() const
{
    return TransformEvaluators<std::string>([](const Evaluator& eval) { return eval.ReadName(); });
}

std::vector<int> Integrator::NumFlipFlops() const
{
    return TransformEvaluators<int>([](const Evaluator& eval) { return eval.NumFlipFlops(); });
}

int Integrator::MaxNumFlipFlops() const { return MaxElement<int>(NumFlipFlops()); }

float Integrator::MaxAlphaPopulated() const
{
    const auto functor = [](const Evaluator& eval) {
        return (eval.IsValid() ? eval.Alpha().UsedEntriesRatio() : NEG_FLOAT_INF);
    };
    auto alphaPopulated = TransformEvaluators<float>(functor);
    return MaxElement<float>(alphaPopulated);
}

float Integrator::MaxBetaPopulated() const
{
    const auto functor = [](const Evaluator& eval) {
        return (eval.IsValid() ? eval.Beta().UsedEntriesRatio() : NEG_FLOAT_INF);
    };
    auto betaPopulated = TransformEvaluators<float>(functor);
    return MaxElement<float>(betaPopulated);
}

double Integrator::AvgZScore() const
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

std::vector<double> Integrator::ZScores() const
{
    return TransformEvaluators<double>([](const Evaluator& eval) { return eval.ZScore(); });
}

std::vector<std::pair<double, double>> Integrator::NormalParameters() const
{
    return TransformEvaluators<std::pair<double, double>>(
        [](const Evaluator& eval) { return eval.NormalParameters(); });
}

void Integrator::MaskIntervals(const size_t radius, const double maxErrRate)
{
    for (auto& eval : evals_)
        if (eval) eval.MaskIntervals(radius, maxErrRate);
}

std::vector<Data::State> Integrator::States() const
{
    return TransformEvaluators<Data::State>([](const Evaluator& eval) { return eval.Status(); });
}

std::vector<StrandType> Integrator::StrandTypes() const
{
    return TransformEvaluators<StrandType>([](const Evaluator& eval) { return eval.Strand(); });
}

const Evaluator& Integrator::GetEvaluator(size_t idx) const { return evals_[idx]; }

const AbstractMatrix& Integrator::Alpha(size_t idx) const { return evals_[idx].Alpha(); }

const AbstractMatrix& Integrator::Beta(size_t idx) const { return evals_[idx].Beta(); }

Mutation Integrator::ReverseComplement(const Mutation& mut) const
{
    size_t newStart = TemplateLength() - mut.End();
    if (mut.IsDeletion())
        return Mutation::Deletion(newStart, mut.Length());
    else if (mut.IsInsertion())
        return Mutation::Insertion(newStart, ::PacBio::Data::ReverseComplement(mut.Bases()));
    // IsSubstitution
    return Mutation::Substitution(newStart, ::PacBio::Data::ReverseComplement(mut.Bases()));
}

State Integrator::AddRead(const PacBio::Data::MappedRead& read)
{
    try {
        return AddRead(GetTemplate(read), read);
    } catch (const TemplateTooSmall& e) {
        return State::TEMPLATE_TOO_SMALL;
    }
}

size_t Integrator::TemplateLength() const { return fwdTpl_.length(); }

char Integrator::operator[](const size_t i) const { return fwdTpl_[i]; }

Integrator::operator std::string() const { return fwdTpl_; }

void Integrator::ApplyMutation(const Mutation& fwdMut)
{
    const Mutation revMut(ReverseComplement(fwdMut));

    std::vector<Mutation> fwdMuts = {fwdMut};
    std::vector<Mutation> revMuts = {revMut};

    fwdTpl_ = ::PacBio::Consensus::ApplyMutations(fwdTpl_, &fwdMuts);
    revTpl_ = ::PacBio::Consensus::ApplyMutations(revTpl_, &revMuts);

    for (auto& eval : evals_) {
        if (eval.Strand() == StrandType::FORWARD)
            eval.ApplyMutation(fwdMut);
        else if (eval.Strand() == StrandType::REVERSE)
            eval.ApplyMutation(revMut);
    }

    assert(fwdTpl_.length() == revTpl_.length());
    assert(fwdTpl_ == ::PacBio::Data::ReverseComplement(revTpl_));
}

void Integrator::ApplyMutations(std::vector<Mutation>* fwdMuts)
{
    std::vector<Mutation> revMuts;

    for (auto it = fwdMuts->crbegin(); it != fwdMuts->crend(); ++it)
        revMuts.emplace_back(ReverseComplement(*it));

    fwdTpl_ = ::PacBio::Consensus::ApplyMutations(fwdTpl_, fwdMuts);
    revTpl_ = ::PacBio::Consensus::ApplyMutations(revTpl_, &revMuts);

    for (auto& eval : evals_) {
        if (eval.Strand() == StrandType::FORWARD)
            eval.ApplyMutations(fwdMuts);
        else if (eval.Strand() == StrandType::REVERSE)
            eval.ApplyMutations(&revMuts);
    }

    assert(fwdTpl_.length() == revTpl_.length());
    assert(fwdTpl_ == ::PacBio::Data::ReverseComplement(revTpl_));
}

std::unique_ptr<AbstractTemplate> Integrator::GetTemplate(const PacBio::Data::MappedRead& read)
{
    const size_t len = read.TemplateEnd - read.TemplateStart;

    if (read.Strand == StrandType::FORWARD) {
        const size_t start = read.TemplateStart;
        const size_t end = read.TemplateEnd;

        return std::unique_ptr<AbstractTemplate>(new Template(fwdTpl_.substr(start, len),
                                                              ModelFactory::Create(read), start,
                                                              end, read.PinStart, read.PinEnd));
    } else if (read.Strand == StrandType::REVERSE) {
        const size_t start = revTpl_.size() - read.TemplateEnd;
        const size_t end = revTpl_.size() - read.TemplateStart;

        return std::unique_ptr<AbstractTemplate>(new Template(revTpl_.substr(start, len),
                                                              ModelFactory::Create(read), start,
                                                              end, read.PinEnd, read.PinStart));
    }

    throw std::invalid_argument("read is unmapped!");
}

}  // namespace Consensus
}  // namespace PacBio
