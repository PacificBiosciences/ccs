
#include <cmath>

#include <pacbio/consensus/Evaluator.h>

#include "EvaluatorImpl.h"

namespace PacBio {
namespace Consensus {

constexpr size_t EXTEND_BUFFER_COLUMNS = 8;

Evaluator::Evaluator(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr,
                     const double scoreDiff)
    : impl_{std::unique_ptr<EvaluatorImpl>(
          new EvaluatorImpl(std::forward<std::unique_ptr<AbstractTemplate>>(tpl), mr, scoreDiff))}
{
    if (impl_.get() == nullptr) throw std::runtime_error("failed to construct EvaluatorImpl");
}

Evaluator::Evaluator(Evaluator&& eval) : impl_{std::move(eval.impl_)} {}
Evaluator::~Evaluator() {}
double Evaluator::LL(const Mutation& mut) { return impl_->LL(mut); }
double Evaluator::LL() const { return impl_->LL(); }
std::tuple<double, double> Evaluator::NormalParameters() const { return impl_->NormalParameters(); }
double Evaluator::ZScore() const { return impl_->ZScore(); }
void Evaluator::ApplyMutation(const Mutation& mut) { impl_->ApplyMutation(mut); }
void Evaluator::ApplyMutations(std::vector<Mutation>* muts) { impl_->ApplyMutations(muts); }
/*
Matrix* Evaluator::AlphaMatrix();
Matrix* Evaluator::BetaMatrix();
*/

}  // namespace Consensus
}  // namespace PacBio
