
#pragma once

#include <memory>
#include <utility>
#include <vector>

#include <pacbio/consensus/Evaluator.h>
#include <pacbio/consensus/Read.h>
#include <pacbio/consensus/Template.h>

#include "Recursor.h"
#include "matrix/ScaledMatrix.h"

namespace PacBio {
namespace Consensus {

class EvaluatorImpl
{
public:
    EvaluatorImpl(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr,
                  double scoreDiff = 12.5);

    std::string ReadName() const;

    double LL(const Mutation& mut);
    double LL() const;

    // TODO: Comments are nice!  Explain what this is about---ZScore calculation?
    std::pair<double, double> NormalParameters() const;

    double ZScore() const;

    bool ApplyMutation(const Mutation& mut);
    bool ApplyMutations(std::vector<Mutation>* muts);

private:
    void Recalculate();

private:
    std::unique_ptr<AbstractRecursor>
        recursor_;  // TODO: does this need to be a pointer?  is it always non-null?
                    // are we making it a UP just so we can do a fwd decl and still have
                    // RAII semantics?
    ScaledMatrix alpha_;
    ScaledMatrix beta_;
    ScaledMatrix extendBuffer_;

    friend class Evaluator;
};

}  // namespace Consensus
}  // namespace PacBio
