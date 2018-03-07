// Author: Lance Hepler

#include <cmath>
#include <memory>
#include <string>

#include <pacbio/consensus/Evaluator.h>
#include <pacbio/exception/InvalidEvaluatorException.h>
#include <pacbio/exception/StateError.h>
#include <pbcopper/logging/Logging.h>

#include "Constants.h"
#include "EvaluatorImpl.h"

using namespace PacBio::Data;
using namespace PacBio::Exception;

namespace PacBio {
namespace Consensus {

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
        impl_ = std::make_unique<EvaluatorImpl>(std::move(tpl), mr, scoreDiff);
        CheckZScore(minZScore, mr.Model);
    } catch (const StateError& e) {
        Status(e.WhatState());
    }
}

Evaluator::Evaluator(Evaluator&& eval) : impl_{std::move(eval.impl_)}, curState_{eval.curState_} {}

Evaluator& Evaluator::operator=(Evaluator&& eval)
{
    if (eval == *this) return *this;
    impl_ = std::move(eval.impl_);
    curState_ = eval.curState_;
    return *this;
}

Evaluator::~Evaluator() = default;

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

    if (!IsValid()) {
        return NEG_DBL_INF;
    }

    // for multi-base mutations, invoke the entire machinery
    else if (mut.EditDistance() > 1) {
        boost::optional<MutatedTemplate> mutTpl = impl_->tpl_->Mutate(mut);
        if (!mutTpl) return NEG_DBL_INF;
        auto newTpl = std::make_unique<MutatedTemplate>(std::move(*mutTpl));
        EvaluatorImpl tmp(std::move(newTpl), impl_->recursor_->read_, impl_->recursor_->scoreDiff_);
        ll = tmp.LL();
    }

    // single-base mutations employ the alpha-beta stitching
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

void Evaluator::MaskIntervals(const size_t radius, const double maxErrRate)
{
    if (IsValid()) return impl_->MaskIntervals(radius, maxErrRate);
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
