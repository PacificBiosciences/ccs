
#pragma once

namespace PacBio {
namespace Consensus {

class Evaluator
{
public:
    Evaluator(std::unique_ptr<AbstractTemplate>&& tpl, MappedRead& mr)
        : tpl_{tpl}
        , mr_{mr}
    {
        if (tpl_.get() == nullptr)
            throw std::invalid_argument("template must be non-null");
    }

    double LL(const Mutation& mut)
    {
        tpl_->Mutate(mut);
        // do JoinAlphaBeta thingys
        tpl_->Reset();
    }

    double LL() const
    {
        return beta_(0, 0);
    }

    std::tuple<double, double> NormalParameters() const
    {
        return tpl_->NormalParameters(mr.TemplateStart, mr.TemplateEnd);
    }

    double ZScore() const
    {
        double mean, var;
        std::tie(mean, var) = NormalParameters();
        return (LL() - mean) / std::sqrt(var);
    }

    void ApplyMutation(const Mutation& mut);
    void ApplyMutations(const std::vector<Mutation>& muts);

    Matrix* AlphaMatrix();
    Matrix* BetaMatrix();

private:

    class EvaluatorImpl;
   
    std::unique_ptr<EvaluatorImpl> impl_;
    // todo store smart ptr here instead
    std::unique_ptr<AbstractTemplate> tpl_;
    MappedRead mr_;
};

} // namespace Consensus
} // namespace PacBio
