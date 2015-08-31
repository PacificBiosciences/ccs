
#include <utility>

#include <pacbio/consensus/Integrator.h>

namespace PacBio {
namespace Consensus {

AddReadResult AbstractIntegrator::AddRead(Evaluator&& eval)
{
    try
    {
        evals_.emplace_back(std::forward<Evaluator>(eval));
    }
    catch (AlphaBetaMismatch& e)
    {
        return ALPHA_BETA_MISMATCH;
    }
    // TODO(lhepler): do we really want other?
    catch (...)
    {
        return OTHER;
    }

    if (!std::isnan(cfg_.MinZScore) &&
        evals_.back().ZScore() < cfg_.MinZScore)
    {
        evals_.pop_back();
        return POOR_ZSCORE;
    }

    return SUCCESS;
}

double AbstractIntegrator::LL(const Mutation& mut)
{
    double ll = 0.0;
    for (auto& eval : evals_)
    {
        ll += eval.LL(mut);
    }
    return ll;
}

double AbstractIntegrator::LL() const
{
    double ll = 0.0;
    for (const auto& eval : evals_)
    {
        ll += eval.LL();
    }
    return ll;
}

double AbstractIntegrator::AvgZScore() const
{
    double mean = 0.0, var = 0.0;
    for (const auto& eval : evals_)
    {
        double m, v;
        std::tie(m, v) = eval.NormalParameters();
        mean += m;
        var  += v;
    }
    const size_t n = evals_.size();
    return (LL()/n - mean/n) / std::sqrt(var/n);
}

MonoMolecularIntegrator::MonoMolecularIntegrator(const std::string& tpl,
                                                 const IntegratorConfig& cfg,
                                                 const SNR& snr,
                                                 const std::string& model)
    : AbstractIntegrator(cfg)
    , mdl_{model}
    , tpl_(tpl, cfg_.ParamTable->At(mdl_, snr))
{ }

AddReadResult MonoMolecularIntegrator::AddRead(const MappedRead& read)
{
    if (read.Model != mdl_)
        return OTHER;
  
    return AbstractIntegrator::AddRead(Evaluator(
                std::unique_ptr<AbstractTemplate>(new VirtualTemplate(tpl_, read.TemplateStart, read.TemplateEnd)),
                read));
}

double MonoMolecularIntegrator::LL(const Mutation& mut)
{
    tpl_.Mutate(mut);
    const double ll = AbstractIntegrator::LL(mut);
    tpl_.Reset();
    return ll;
}

MultiMolecularIntegrator::MultiMolecularIntegrator(const std::string& tpl,
                                                   const IntegratorConfig& cfg)
    : AbstractIntegrator(cfg)
    , tpl_{tpl}
{ }

AddReadResult MultiMolecularIntegrator::AddRead(const MappedRead& read, const SNR& snr)
{
    return AbstractIntegrator::AddRead(Evaluator(
                std::unique_ptr<AbstractTemplate>(
                    new Template(tpl_, cfg_.ParamTable->At(read.Model, snr))),
                read));
}

} // namespace Consensus
} // namespace PacBio
