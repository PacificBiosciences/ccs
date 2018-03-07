// Author: Brett Bowman

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

template <bool AllowInvalidEvaluators>
inline double Integrator::SingleEvaluatorLL(Evaluator* const eval, const Mutation& fwdMut) const
{
    const Mutation revMut{ReverseComplement(fwdMut)};
    double result{AllowInvalidEvaluators ? NEG_DBL_INF : 0};

    switch (eval->Strand()) {
        case StrandType::FORWARD:
            result = eval->LL(fwdMut);
            break;
        case StrandType::REVERSE:
            result = eval->LL(revMut);
            break;
        case StrandType::UNMAPPED:
            if (AllowInvalidEvaluators == false)
                throw InvalidEvaluatorException("Unmapped read in mutation testing");
            break;
        default:
            throw std::runtime_error("Unknown StrandType");
    }

    return result;
}

double Integrator::LL(const Mutation& fwdMut)
{
    double ll = 0.0;
    for (auto& e : evals_) {
        // Skip invalid Evaluators
        if (!e.IsValid()) continue;

        ll += SingleEvaluatorLL<false>(&e, fwdMut);
    }
    return ll;
}

double Integrator::LL() const
{
    const auto functor = [](const Evaluator& eval) { return eval.IsValid() ? eval.LL() : 0; };
    const auto lls = TransformEvaluators<double>(functor);
    return std::accumulate(lls.cbegin(), lls.cend(), 0.0);
}

std::vector<double> Integrator::LLs(const Mutation& fwdMut)
{
    // Compute individual LLs of each Evaluator
    std::vector<double> lls;
    lls.reserve(evals_.size());
    for (auto& e : evals_) {
        const double ll = SingleEvaluatorLL<true>(&e, fwdMut);
        lls.emplace_back(ll);
    }
    return lls;
}

std::vector<double> Integrator::LLs() const
{
    const auto functor = [](const Evaluator& eval) { return eval.LL(); };
    return TransformEvaluators<double>(functor);
}

std::array<std::pair<char, int>, 4> Integrator::BestMutationHistogram(const size_t start,
                                                                      const MutationType mutType)
{
    // Histogram makes no sense for deletions
    if (mutType == MutationType::DELETION)
        throw std::runtime_error("Cannot create a histogram over a deletion mutation");

    std::array<std::pair<char, int>, 4> result{{{'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}}};

    for (auto& eval : evals_) {
        // Skip invalid Evaluators
        if (!eval.IsValid()) continue;

        double curBestLL = eval.LL();
        int8_t curBestIdx = -1;  // -1 is a sentinel for the current template being the best so far

        for (int8_t i = 0; i < 4; ++i) {
            const char c = result[i].first;
            const Mutation testMut{(mutType == MutationType::SUBSTITUTION)
                                       ? Mutation::Substitution(start, c)
                                       : Mutation::Insertion(start, c)};

            const double ll = SingleEvaluatorLL<false>(&eval, testMut);

            if (curBestLL < ll) {
                curBestLL = ll;
                curBestIdx = i;
            }
        }

        if (curBestIdx != -1) ++result[curBestIdx].second;
    }

    std::sort(result.begin(), result.end(),
              [](const std::pair<char, int>& lhs, const std::pair<char, int>& rhs) {
                  return rhs.second < lhs.second;
              });

    return result;
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

        return std::make_unique<Template>(fwdTpl_.substr(start, len), ModelFactory::Create(read),
                                          start, end, read.PinStart, read.PinEnd);
    } else if (read.Strand == StrandType::REVERSE) {
        const size_t start = revTpl_.size() - read.TemplateEnd;
        const size_t end = revTpl_.size() - read.TemplateStart;

        return std::make_unique<Template>(revTpl_.substr(start, len), ModelFactory::Create(read),
                                          start, end, read.PinEnd, read.PinStart);
    }

    throw std::invalid_argument("read is unmapped!");
}

}  // namespace Consensus
}  // namespace PacBio
