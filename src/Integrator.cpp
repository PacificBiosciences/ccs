
#include <cmath>
#include <utility>

#include <pacbio/consensus/Integrator.h>
#include <pacbio/consensus/Sequence.h>

#include "ModelFactory.h"

namespace PacBio {
namespace Consensus {

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
AddReadResult AbstractIntegrator::AddRead(std::unique_ptr<AbstractTemplate>&& tpl,
                                          const MappedRead& read)
{
    std::unique_ptr<Evaluator> eval(nullptr);
    AddReadResult result = AddReadResult::SUCCESS;

    try {
        eval.reset(new Evaluator(std::move(tpl), read, cfg_.ScoreDiff));
        const double zScore = eval->ZScore();
        if (!std::isnan(cfg_.MinZScore) && (!std::isfinite(zScore) || zScore < cfg_.MinZScore)) {
            eval.reset(nullptr);
            result = AddReadResult::POOR_ZSCORE;
        }
    } catch (AlphaBetaMismatch& e) {
        eval.reset(nullptr);
        result = AddReadResult::ALPHA_BETA_MISMATCH;
    }
    // TODO(lhepler): do we really want other?
    catch (...) {
        eval.reset(nullptr);
        result = AddReadResult::OTHER;
    }

    evals_.emplace_back(ReadState(std::move(eval), read.Strand));

    return result;
}

double AbstractIntegrator::LL(const Mutation& fwdMut)
{
    const Mutation revMut(ReverseComplement(fwdMut));
    double ll = 0.0;
    for (auto& eval : evals_) {
        if (eval) {
            if (eval.Strand == StrandEnum::FORWARD)
                ll += eval->LL(fwdMut);
            else
                ll += eval->LL(revMut);
        }
    }
    return ll;
}

double AbstractIntegrator::LL() const
{
    double ll = 0.0;
    for (const auto& eval : evals_) {
        if (eval) ll += eval->LL();
    }
    return ll;
}

double AbstractIntegrator::AvgZScore() const
{
    double mean = 0.0, var = 0.0;
    size_t n = 0;
    for (const auto& eval : evals_) {
        if (eval) {
            double m, v;
            std::tie(m, v) = eval->NormalParameters();
            mean += m;
            var += v;
            ++n;
        }
    }
    return (LL() / n - mean / n) / std::sqrt(var / n);
}

std::vector<double> AbstractIntegrator::ZScores() const
{
    std::vector<double> results;
    results.reserve(evals_.size());
    for (const auto& eval : evals_) {
        if (eval) {
            double mean, var;
            std::tie(mean, var) = eval->NormalParameters();
            results.emplace_back((eval->LL() - mean) / std::sqrt(var));
        } else
            results.emplace_back(std::numeric_limits<double>::quiet_NaN());
    }
    return results;
}

Mutation AbstractIntegrator::ReverseComplement(const Mutation& mut) const
{
    return Mutation(mut.Type, Length() - mut.End(), Complement(mut.Base));
}

AbstractIntegrator::ReadState::ReadState(std::unique_ptr<Evaluator>&& ptr, StrandEnum strand)
    : std::unique_ptr<Evaluator>(std::forward<std::unique_ptr<Evaluator>>(ptr)), Strand(strand)
{
}

MonoMolecularIntegrator::MonoMolecularIntegrator(const std::string& tpl,
                                                 const IntegratorConfig& cfg, const SNR& snr,
                                                 const std::string& model)
    : AbstractIntegrator(cfg)
    , mdl_{model}
    , fwdTpl_(tpl, ModelFactory::Create(mdl_, snr))
    , revTpl_(::PacBio::Consensus::ReverseComplement(tpl), ModelFactory::Create(mdl_, snr))
{
}

MonoMolecularIntegrator::MonoMolecularIntegrator(MonoMolecularIntegrator&& mmi)
    : AbstractIntegrator(std::move(mmi))
    , mdl_{mmi.mdl_}
    , fwdTpl_{std::move(mmi.fwdTpl_)}
    , revTpl_{std::move(mmi.revTpl_)}
{
}

AddReadResult MonoMolecularIntegrator::AddRead(const MappedRead& read)
{
    if (read.Model != mdl_) return AddReadResult::OTHER;

    if (read.Strand == StrandEnum::FORWARD)
        return AbstractIntegrator::AddRead(
            std::unique_ptr<AbstractTemplate>(new VirtualTemplate(
                fwdTpl_, read.TemplateStart, read.TemplateEnd, read.PinStart, read.PinEnd)),
            read);

    return AbstractIntegrator::AddRead(
        std::unique_ptr<AbstractTemplate>(new VirtualTemplate(revTpl_, Length() - read.TemplateEnd,
                                                              Length() - read.TemplateStart,
                                                              read.PinEnd, read.PinStart)),
        read);
}

size_t MonoMolecularIntegrator::Length() const { return fwdTpl_.TrueLength(); }
char MonoMolecularIntegrator::operator[](const size_t i) const { return fwdTpl_[i].Base; }
MonoMolecularIntegrator::operator std::string() const
{
    std::string result;

    result.resize(fwdTpl_.Length());

    for (size_t i = 0; i < fwdTpl_.Length(); ++i)
        result[i] = fwdTpl_[i].Base;

    return result;
}

double MonoMolecularIntegrator::LL(const Mutation& fwdMut)
{
    const Mutation revMut(ReverseComplement(fwdMut));
    fwdTpl_.Mutate(fwdMut);
    revTpl_.Mutate(revMut);
    const double ll = AbstractIntegrator::LL(fwdMut);
    fwdTpl_.Reset();
    revTpl_.Reset();
    return ll;
}

void MonoMolecularIntegrator::ApplyMutation(const Mutation& fwdMut)
{
    const Mutation revMut(ReverseComplement(fwdMut));

    fwdTpl_.ApplyMutation(fwdMut);
    revTpl_.ApplyMutation(revMut);

    for (auto& eval : evals_) {
        if (eval) {
            if (eval.Strand == StrandEnum::FORWARD)
                eval->ApplyMutation(fwdMut);
            else
                eval->ApplyMutation(revMut);
        }
    }

    assert(fwdTpl_.Length() == revTpl_.Length());

#ifndef NDEBUG
    std::string fwd;
    std::string rev;

    for (size_t i = 0; i < Length(); ++i) {
        fwd.push_back(fwdTpl_[i].Base);
        rev.push_back(revTpl_[i].Base);
    }

#endif

    assert(fwd == ::PacBio::Consensus::ReverseComplement(rev));
}

void MonoMolecularIntegrator::ApplyMutations(std::vector<Mutation>* fwdMuts)
{
    std::vector<Mutation> revMuts;

    for (auto it = fwdMuts->crbegin(); it != fwdMuts->crend(); ++it)
        revMuts.emplace_back(ReverseComplement(*it));

    fwdTpl_.ApplyMutations(fwdMuts);
    revTpl_.ApplyMutations(&revMuts);

    for (auto& eval : evals_) {
        if (eval) {
            if (eval.Strand == StrandEnum::FORWARD)
                eval->ApplyMutations(fwdMuts);
            else
                eval->ApplyMutations(&revMuts);
        }
    }

    assert(fwdTpl_.Length() == revTpl_.Length());

#ifndef NDEBUG
    std::string fwd;
    std::string rev;

    for (size_t i = 0; i < Length(); ++i) {
        fwd.push_back(fwdTpl_[i].Base);
        rev.push_back(revTpl_[i].Base);
    }
#endif

    assert(fwd == ::PacBio::Consensus::ReverseComplement(rev));
}

MultiMolecularIntegrator::MultiMolecularIntegrator(const std::string& tpl,
                                                   const IntegratorConfig& cfg)
    : AbstractIntegrator(cfg), fwdTpl_{tpl}, revTpl_{::PacBio::Consensus::ReverseComplement(tpl)}
{
}

AddReadResult MultiMolecularIntegrator::AddRead(const MappedRead& read, const SNR& snr)
{
    const size_t len = read.TemplateEnd - read.TemplateStart;

    if (read.Strand == StrandEnum::FORWARD) {
        const size_t start = read.TemplateStart;
        const size_t end = read.TemplateEnd;

        return AbstractIntegrator::AddRead(
            std::unique_ptr<AbstractTemplate>(new Template(fwdTpl_.substr(start, len),
                                                           ModelFactory::Create(read.Model, snr),
                                                           start, end, read.PinStart, read.PinEnd)),
            read);
    } else {
        const size_t start = revTpl_.size() - read.TemplateEnd;
        const size_t end = revTpl_.size() - read.TemplateStart;

        return AbstractIntegrator::AddRead(
            std::unique_ptr<AbstractTemplate>(new Template(revTpl_.substr(start, len),
                                                           ModelFactory::Create(read.Model, snr),
                                                           start, end, read.PinEnd, read.PinStart)),
            read);
    }

    return AddReadResult::OTHER;
}

size_t MultiMolecularIntegrator::Length() const { return fwdTpl_.length(); }
char MultiMolecularIntegrator::operator[](const size_t i) const { return fwdTpl_[i]; }
MultiMolecularIntegrator::operator std::string() const { return fwdTpl_; }
void MultiMolecularIntegrator::ApplyMutation(const Mutation& fwdMut)
{
    const Mutation revMut(ReverseComplement(fwdMut));

    std::vector<Mutation> fwdMuts = {fwdMut};
    std::vector<Mutation> revMuts = {revMut};

    fwdTpl_ = ::PacBio::Consensus::ApplyMutations(fwdTpl_, &fwdMuts);
    revTpl_ = ::PacBio::Consensus::ApplyMutations(revTpl_, &revMuts);

    for (auto& eval : evals_) {
        if (eval) {
            if (eval.Strand == StrandEnum::FORWARD)
                eval->ApplyMutation(fwdMut);
            else
                eval->ApplyMutation(revMut);
        }
    }

    assert(fwdTpl_.length() == revTpl_.length());
    assert(fwdTpl_ == ::PacBio::Consensus::ReverseComplement(revTpl_));
}

void MultiMolecularIntegrator::ApplyMutations(std::vector<Mutation>* fwdMuts)
{
    std::vector<Mutation> revMuts;

    for (auto it = fwdMuts->crbegin(); it != fwdMuts->crend(); ++it)
        revMuts.emplace_back(ReverseComplement(*it));

    fwdTpl_ = ::PacBio::Consensus::ApplyMutations(fwdTpl_, fwdMuts);
    revTpl_ = ::PacBio::Consensus::ApplyMutations(revTpl_, &revMuts);

    for (auto& eval : evals_) {
        if (eval) {
            if (eval.Strand == StrandEnum::FORWARD)
                eval->ApplyMutations(fwdMuts);
            else
                eval->ApplyMutations(&revMuts);
        }
    }

    assert(fwdTpl_.length() == revTpl_.length());
    assert(fwdTpl_ == ::PacBio::Consensus::ReverseComplement(revTpl_));
}

}  // namespace Consensus
}  // namespace PacBio
