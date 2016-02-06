
#pragma once

#include <memory>
#include <utility>
#include <vector>

#include <pacbio/consensus/Read.h>
#include <pacbio/consensus/Template.h>

#include "matrix/ScaledMatrix.h"
#include "Recursor.h"

namespace PacBio {
namespace Consensus {

class EvaluatorImpl
{
public:
    EvaluatorImpl(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr,
                  double scoreDiff = 12.5);

    double LL(const Mutation& mut);
    double LL() const;

    std::pair<double, double> NormalParameters() const;

    double ZScore() const;

    bool ApplyMutation(const Mutation& mut);
    bool ApplyMutations(std::vector<Mutation>* muts);

private:
    void Recalculate();

private:
    Recursor recursor_;

    ScaledMatrix alpha_;
    ScaledMatrix beta_;
    ScaledMatrix extendBuffer_;
};

}  // namespace Consensus
}  // namespace PacBio
