// Author: Lance Hepler

#pragma once

#include <memory>
#include <utility>
#include <vector>

#include <pacbio/consensus/Evaluator.h>
#include <pacbio/consensus/IntervalMask.h>
#include <pacbio/consensus/MatrixViewConvention.h>
#include <pacbio/consensus/Template.h>
#include <pacbio/data/Read.h>

#include "Recursor.h"
#include "matrix/ScaledMatrix.h"

namespace PacBio {
namespace Consensus {

class EvaluatorImpl
{
public:
    EvaluatorImpl(std::unique_ptr<AbstractTemplate>&& tpl, const PacBio::Data::MappedRead& mr,
                  double scoreDiff = 12.5);

    std::string ReadName() const;

    double LL(const Mutation& mut);
    double LL() const;

    // Interval masking methods
    void MaskIntervals(size_t radius, double maxErrRate);

    // TODO: Comments are nice!  Explain what this is about---ZScore calculation?
    std::pair<double, double> NormalParameters() const;

    double ZScore() const;

    bool ApplyMutation(const Mutation& mut);
    bool ApplyMutations(std::vector<Mutation>* muts);

    int NumFlipFlops() const { return numFlipFlops_; }

public:
    const AbstractMatrix& Alpha() const;
    const AbstractMatrix& Beta() const;

public:
    const AbstractMatrix* AlphaView(MatrixViewConvention c) const;
    const AbstractMatrix* BetaView(MatrixViewConvention c) const;

private:
    void Recalculate();

private:
    std::unique_ptr<AbstractTemplate> tpl_;
    std::unique_ptr<AbstractRecursor> recursor_;
    IntervalMask mask_;
    ScaledMatrix alpha_;
    ScaledMatrix beta_;
    ScaledMatrix extendBuffer_;

    int numFlipFlops_;

    friend class Evaluator;
};

}  // namespace Consensus
}  // namespace PacBio
