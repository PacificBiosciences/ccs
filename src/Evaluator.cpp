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
#include <pacbio/exception/InvalidEvaluatorException.h>
#include <pacbio/exception/StateError.h>
#include <pbcopper/logging/Logging.h>

#include "EvaluatorImpl.h"

using namespace PacBio::Data;
using namespace PacBio::Exception;

namespace PacBio {
namespace Consensus {
namespace {

constexpr double NEG_DBL_INF = -std::numeric_limits<double>::infinity();
constexpr int NEG_INT_INF = -std::numeric_limits<int>::infinity();
constexpr float NEG_FLOAT_INF = -std::numeric_limits<float>::infinity();
constexpr size_t EXTEND_BUFFER_COLUMNS = 8;

}  // anonymous namespace

Evaluator::Evaluator(const State state) : impl_{nullptr}, curState_{state}
{
    if (curState_ == State::VALID)
        throw std::invalid_argument("cannot initialize a dummy Evaluator with VALID state");
}

Evaluator::Evaluator(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr,
                     const double minZScore, const double scoreDiff)
    : impl_{nullptr}, curState_{State::VALID}
{
    try {
        impl_.reset(
            new EvaluatorImpl(std::forward<std::unique_ptr<AbstractTemplate>>(tpl), mr, scoreDiff));
        CheckZScore(minZScore, mr.Model);
    } catch (const StateError& e) {
        Status(e.WhatState());
    }
}

Evaluator::Evaluator(Evaluator&& eval) : impl_{std::move(eval.impl_)}, curState_{eval.curState_} {}

Evaluator& Evaluator::operator=(Evaluator&& eval)
{
    if (eval == *this) return *this;
    impl_ = move(eval.impl_);
    curState_ = eval.curState_;
    return *this;
}

Evaluator::~Evaluator() {}

size_t Evaluator::Length() const
{
    if (IsValid()) return impl_->tpl_->Length();
    return 0;
}

StrandType Evaluator::Strand() const
{
    if (IsValid()) return impl_->recursor_->read_.Strand;
    return StrandType::UNMAPPED;
}

std::string Evaluator::ReadName() const
{
    if (IsValid()) return impl_->ReadName();
    return "*Inactive evaluator*";
}

double Evaluator::LL(const Mutation& mut)
{
    double ll;

    if (!IsValid())
        return NEG_DBL_INF;

    else if (mut.EditDistance() > 1) {
        boost::optional<MutatedTemplate> mutTpl = impl_->tpl_->Mutate(mut);
        if (!mutTpl) return NEG_DBL_INF;
        std::unique_ptr<AbstractTemplate> newTpl(new MutatedTemplate(std::move(*mutTpl)));
        EvaluatorImpl tmp(std::move(newTpl), impl_->recursor_->read_, impl_->recursor_->scoreDiff_);
        ll = tmp.LL();
    }

    // single-base mutation
    else {
        ll = impl_->LL(mut);
    }

    // If the mutation of interest caused a corner-case failure,
    // release this Evaluator and report this issue via an exception.
    if (std::isinf(ll)) {
        const std::string name = ReadName();
        Invalidate();
        throw InvalidEvaluatorException("negative inf in mutation testing: '" + name + "'");
    }

    return ll;
}

double Evaluator::LL() const
{
    if (IsValid()) return impl_->LL();
    return NEG_DBL_INF;
}

std::pair<double, double> Evaluator::NormalParameters() const
{
    if (IsValid()) return impl_->NormalParameters();
    return std::make_pair(NEG_DBL_INF, NEG_DBL_INF);
}

double Evaluator::ZScore() const
{
    if (IsValid()) return impl_->ZScore();
    return NEG_DBL_INF;
}

int Evaluator::NumFlipFlops() const
{
    if (IsValid()) return impl_->NumFlipFlops();
    return NEG_INT_INF;
}

bool Evaluator::ApplyMutation(const Mutation& mut)
{
    bool mutApplied = false;
    if (IsValid()) {
        try {
            mutApplied = impl_->ApplyMutation(mut);
        } catch (const StateError& e) {
            Status(e.WhatState());
        }
    }
    return mutApplied;
}

bool Evaluator::ApplyMutations(std::vector<Mutation>* muts)
{
    bool mutsApplied = false;
    if (IsValid()) {
        try {
            mutsApplied = impl_->ApplyMutations(muts);
        } catch (const StateError& e) {
            Status(e.WhatState());
        }
    }
    return mutsApplied;
}

void Evaluator::Status(State nextState)
{
    // Allow transition from VALID to anything and
    // from anything to DISABLED
    if (curState_ == State::VALID)
        curState_ = nextState;
    else
        PBLOG_ERROR << "Log this behaviour and return";

    if (curState_ != State::VALID) impl_.reset(nullptr);
}

void Evaluator::Release() { Status(State::MANUALLY_RELEASED); }

void Evaluator::Invalidate() { Status(State::INVALID); }

const AbstractMatrix& Evaluator::Alpha() const
{
    if (IsValid()) {
        return impl_->Alpha();
    } else {
        return ScaledMatrix::Null();
    }
}

const AbstractMatrix& Evaluator::Beta() const
{
    if (IsValid()) {
        return impl_->Beta();
    } else {
        return ScaledMatrix::Null();
    }
}

const AbstractMatrix* Evaluator::AlphaView(MatrixViewConvention c) const
{
    if (IsValid()) {
        return impl_->AlphaView(c);
    } else {
        return nullptr;
    }
}

const AbstractMatrix* Evaluator::BetaView(MatrixViewConvention c) const
{
    if (IsValid()) {
        return impl_->BetaView(c);
    } else {
        return nullptr;
    }
}

void Evaluator::CheckZScore(const double minZScore, const std::string& model)
{
    // the zscore filter is disabled under the following conditions
    // - unsupported model (anything not P6-C4)
    if (model.find("P6-C4") == std::string::npos) return;

    // - threshold undefined or too low
    if (std::isnan(minZScore) || minZScore <= -100.0) return;

    const double zScore = impl_->ZScore();
    // TODO(lhepler): re-enable this check when the zscore bits are working again
    // assert(std::isfinite(zScore));
    if (!std::isfinite(zScore) || zScore < minZScore) Status(State::POOR_ZSCORE);
}
}  // namespace Consensus
}  // namespace PacBio
