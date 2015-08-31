
#pragma once

#include <cmath>
#include <memory>
#include <vector>

#include <pacbio/consensus/Read.h>
#include <pacbio/consensus/Template.h>

namespace PacBio {
namespace Consensus {

class Evaluator
{
public:
    Evaluator(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr)
        : tpl_(std::forward<std::unique_ptr<AbstractTemplate>>(tpl))
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
        return 0.0;
    }

    double LL() const
    {
        return 0.0;
    }

    std::tuple<double, double> NormalParameters() const
    {
        return tpl_->NormalParameters(mr_.TemplateStart, mr_.TemplateEnd);
    }

    double ZScore() const
    {
        double mean, var;
        std::tie(mean, var) = NormalParameters();
        return (LL() - mean) / std::sqrt(var);
    }

    void ApplyMutation(const Mutation& mut);
    void ApplyMutations(const std::vector<Mutation>& muts);

    /*
    Matrix* AlphaMatrix();
    Matrix* BetaMatrix();
    */

private:
    std::unique_ptr<AbstractTemplate> tpl_;
    MappedRead mr_;
};

} // namespace Consensus
} // namespace PacBio
