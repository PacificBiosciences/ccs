// Copyright (c) 2011-2016, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

#include <cmath>
#include <memory>
#include <string>

#include <pacbio/consensus/Evaluator.h>
#include <pacbio/consensus/Exceptions.h>

#include "EvaluatorImpl.h"

namespace PacBio {
namespace Consensus {
namespace {

constexpr double NEG_INF = -std::numeric_limits<double>::infinity();
constexpr size_t EXTEND_BUFFER_COLUMNS = 8;

}  // anonymous namespace

Evaluator::Evaluator(const EvaluatorState state) : impl_{nullptr}, state_{state}
{
    if (state_ == EvaluatorState::VALID)
        throw std::invalid_argument("cannot initialize a dummy Evaluator with VALID state");
}

Evaluator::Evaluator(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr,
                     const double minZScore, const double scoreDiff)
    : impl_{nullptr}, state_{EvaluatorState::VALID}
{
    try {
        impl_.reset(
            new EvaluatorImpl(std::forward<std::unique_ptr<AbstractTemplate>>(tpl), mr, scoreDiff));

        const double zScore = impl_->ZScore();

        // the zscore filter is disabled under the following conditions
        if ((mr.Model.find("S/P1-C1") != std::string::npos) ||
            (mr.Model.find("S/P2-C2/prospective-compatible") != std::string::npos))
            goto end;
        if (minZScore <= -100.0) goto end;
        if (std::isnan(minZScore)) goto end;

        // TODO(lhepler): re-enable this check when the zscore bits are working again
        // assert(std::isfinite(zScore));
        if (!std::isfinite(zScore) || zScore < minZScore) state_ = EvaluatorState::POOR_ZSCORE;
    } catch (AlphaBetaMismatch&) {
        state_ = EvaluatorState::ALPHA_BETA_MISMATCH;
    }

end:
    CheckInvariants();
}

Evaluator::Evaluator(Evaluator&& eval) : impl_{std::move(eval.impl_)}, state_{eval.state_}
{
    CheckInvariants();
}

Evaluator& Evaluator::operator=(Evaluator&& eval)
{
    impl_ = move(eval.impl_);
    state_ = eval.state_;
    return *this;
}

Evaluator::~Evaluator() {}
size_t Evaluator::Length() const
{
    if (impl_) return impl_->recursor_->tpl_->Length();
    return 0;
}

StrandEnum Evaluator::Strand() const
{
    if (impl_) return impl_->recursor_->read_.Strand;
    return StrandEnum::UNMAPPED;
}

Evaluator::operator bool() const { return state_ == EvaluatorState::VALID; }
std::string Evaluator::ReadName() const
{
    if (impl_) {
        return impl_->ReadName();
    } else {
        return "*Inactive evaluator*";
    }
}

double Evaluator::LL(const Mutation& mut)
{
    if (impl_) return impl_->LL(mut);
    return NEG_INF;
}

double Evaluator::LL() const
{
    if (impl_) return impl_->LL();
    return NEG_INF;
}

std::pair<double, double> Evaluator::NormalParameters() const
{
    if (impl_) return impl_->NormalParameters();
    return std::make_pair(NEG_INF, NEG_INF);
}

double Evaluator::ZScore() const
{
    if (impl_) return impl_->ZScore();
    return NEG_INF;
}

bool Evaluator::ApplyMutation(const Mutation& mut)
{
    if (!impl_) return false;
    const bool mutApplied = impl_->ApplyMutation(mut);
    CheckInvariants();
    return mutApplied;
}

bool Evaluator::ApplyMutations(std::vector<Mutation>* muts)
{
    if (!impl_) return false;
    const bool mutsApplied = impl_->ApplyMutations(muts);
    CheckInvariants();
    return mutsApplied;
}

EvaluatorState Evaluator::Status() const { return state_; }
void Evaluator::Release()
{
    state_ = EvaluatorState::DISABLED;
    impl_.reset(nullptr);
}

// TODO: this is no longer a simple "invariants check" function---a CheckInvariants function
// should be const and only do asserts on whether the state is valid.  Rename or refactor this.  The
// state in this class
// is a bit too complex for my taste.
void Evaluator::CheckInvariants()
{
    if (!impl_) return;
    if (impl_->recursor_->tpl_->Length() < 2) state_ = EvaluatorState::NULL_TEMPLATE;
    if (state_ != EvaluatorState::VALID) impl_.reset(nullptr);
}

}  // namespace Consensus
}  // namespace PacBio
