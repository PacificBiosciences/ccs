
#include <cmath>
#include <utility>

#include <pacbio/consensus/Integrator.h>

#include "ModelFactory.h"

namespace PacBio {
namespace Consensus {

IntegratorConfig::IntegratorConfig(const double minZScore, const double scoreDiff)
    : MinZScore{minZScore}, ScoreDiff{scoreDiff}
{
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
    std::unique_ptr<Evaluator> eval;
    AddReadResult result = AddReadResult::SUCCESS;

    try {
        eval.reset(new Evaluator(std::move(tpl), read, cfg_.ScoreDiff));
    } catch (AlphaBetaMismatch& e) {
        eval.reset(nullptr);
        result = AddReadResult::ALPHA_BETA_MISMATCH;
    }
    // TODO(lhepler): do we really want other?
    catch (...) {
        eval.reset(nullptr);
        result = AddReadResult::OTHER;
    }

    if (!std::isnan(cfg_.MinZScore) && eval && eval->ZScore() < cfg_.MinZScore) {
        eval.reset(nullptr);
        result = AddReadResult::POOR_ZSCORE;
    }

    evals_.emplace_back(std::move(eval));

    return result;
}

double AbstractIntegrator::LL(const Mutation& mut)
{
    double ll = 0.0;
    for (auto& eval : evals_) {
        if (eval) ll += eval->LL(mut);
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
    for (const auto& eval : evals_) {
        if (eval) {
            double m, v;
            std::tie(m, v) = eval->NormalParameters();
            mean += m;
            var += v;
        }
    }
    const size_t n = evals_.size();
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

MonoMolecularIntegrator::MonoMolecularIntegrator(const std::string& tpl,
                                                 const IntegratorConfig& cfg, const SNR& snr,
                                                 const std::string& model)
    : AbstractIntegrator(cfg), mdl_{model}, tpl_(tpl, ModelFactory::Create(mdl_, snr))
{
}

MonoMolecularIntegrator::MonoMolecularIntegrator(MonoMolecularIntegrator&& mmi)
    : AbstractIntegrator(std::move(mmi)), mdl_{mmi.mdl_}, tpl_{std::move(mmi.tpl_)}
{
}

AddReadResult MonoMolecularIntegrator::AddRead(const MappedRead& read)
{
    if (read.Model != mdl_) return AddReadResult::OTHER;

    return AbstractIntegrator::AddRead(
        std::unique_ptr<AbstractTemplate>(new VirtualTemplate(
            tpl_, read.TemplateStart, read.TemplateEnd, read.PinStart, read.PinEnd)),
        read);
}

size_t MonoMolecularIntegrator::Length() const { return tpl_.Length(); }
char MonoMolecularIntegrator::operator[](const size_t i) const { return tpl_[i].Base; }
MonoMolecularIntegrator::operator std::string() const
{
    std::string result;

    result.resize(tpl_.Length());

    for (size_t i = 0; i < tpl_.Length(); ++i) {
        result[i] = tpl_[i].Base;
    }

    return result;
}

double MonoMolecularIntegrator::LL(const Mutation& mut)
{
    tpl_.Mutate(mut);
    const double ll = AbstractIntegrator::LL(mut);
    tpl_.Reset();
    return ll;
}

void MonoMolecularIntegrator::ApplyMutation(const Mutation& mut)
{
    tpl_.ApplyMutation(mut);
    for (auto& eval : evals_)
        if (eval) eval->ApplyMutation(mut);
}

void MonoMolecularIntegrator::ApplyMutations(std::vector<Mutation>* muts)
{
    tpl_.ApplyMutations(muts);
    for (auto& eval : evals_)
        if (eval) eval->ApplyMutations(muts);
}

MultiMolecularIntegrator::MultiMolecularIntegrator(const std::string& tpl,
                                                   const IntegratorConfig& cfg)
    : AbstractIntegrator(cfg), tpl_{tpl}
{
}

AddReadResult MultiMolecularIntegrator::AddRead(const MappedRead& read, const SNR& snr)
{
    return AbstractIntegrator::AddRead(
        std::unique_ptr<AbstractTemplate>(new Template(tpl_, ModelFactory::Create(read.Model, snr),
                                                       read.TemplateStart, read.TemplateEnd,
                                                       read.PinStart, read.PinEnd)),
        read);
}

size_t MultiMolecularIntegrator::Length() const { return tpl_.length(); }
char MultiMolecularIntegrator::operator[](const size_t i) const { return tpl_[i]; }
MultiMolecularIntegrator::operator std::string() const { return tpl_; }
void MultiMolecularIntegrator::ApplyMutation(const Mutation& mut)
{
    std::vector<Mutation> muts(1, mut);
    tpl_ = ::PacBio::Consensus::ApplyMutations(tpl_, &muts);
    for (auto& eval : evals_)
        if (eval) eval->ApplyMutation(mut);
}

void MultiMolecularIntegrator::ApplyMutations(std::vector<Mutation>* muts)
{
    tpl_ = ::PacBio::Consensus::ApplyMutations(tpl_, muts);
    for (auto& eval : evals_)
        if (eval) eval->ApplyMutations(muts);
}

}  // namespace Consensus
}  // namespace PacBio
