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

#include <boost/optional.hpp>

#include "EvaluatorImpl.h"

namespace PacBio {
namespace Consensus {
namespace {  // anonymous

constexpr size_t EXTEND_BUFFER_COLUMNS = 8;

#if 0
std::ostream& operator<<(std::ostream& out, const std::pair<size_t, size_t>& x)
{
    return out << '(' << x.first << ", " << x.second << ')';
}

void WriteMatrix(const ScaledMatrix& mat)
{
    std::cerr << std::pair<size_t, size_t>(mat.Rows(), mat.Columns()) << std::endl;

    for (size_t j = 0; j < mat.Columns(); ++j)
        std::cerr << " " << mat.UsedRowRange(j);
    std::cerr << std::endl;

    std::cerr << "lg: ";
    for (size_t j = 0; j < mat.Columns(); ++j)
        std::cerr << "\t" << std::fixed << std::setprecision(3) << mat.GetLogScale(j);
    std::cerr << std::endl;

    std::cerr << "lgS: ";
    double lgS = 0.0;
    for (size_t j = 0; j < mat.Columns(); ++j)
        std::cerr << "\t" << std::fixed << std::setprecision(3) << (lgS += mat.GetLogScale(j));
    std::cerr << std::endl;

    for (size_t i = 0; i < mat.Rows(); ++i)
    {
        for (size_t j = 0; j < mat.Columns(); ++j)
        {
            std::cerr << "\t" << std::fixed << std::setprecision(3) << std::log(mat.Get(i, j)) + mat.GetLogScale(j);
        }
        std::cerr << std::endl;
    }
}
#endif

}  // namespace anonymous

EvaluatorImpl::EvaluatorImpl(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr,
                             const double scoreDiff)
    : recursor_{tpl->CreateRecursor(std::move(tpl), mr, scoreDiff)}
    , alpha_(mr.Length() + 1, recursor_->tpl_->Length() + 1, ScaledMatrix::FORWARD)
    , beta_(mr.Length() + 1, recursor_->tpl_->Length() + 1, ScaledMatrix::REVERSE)
    , extendBuffer_(mr.Length() + 1, EXTEND_BUFFER_COLUMNS, ScaledMatrix::FORWARD)
{
    recursor_->FillAlphaBeta(alpha_, beta_);
}

std::string EvaluatorImpl::ReadName() const { return recursor_->read_.Name; }
double EvaluatorImpl::LL(const Mutation& mut_)
{
    // apply the virtual mutation
    boost::optional<Mutation> mut(recursor_->tpl_->Mutate(mut_));

    // if the resulting template is 0, simulate NULL_TEMPLATE (removal)
    if (recursor_->tpl_->Length() == 0) {
        recursor_->tpl_->Reset();
        return 0.0;
    }

    // if the mutation didn't hit this read, just return the ll as-is
    if (!mut) return LL();

    size_t betaLinkCol = 1 + mut->End();
    size_t absoluteLinkColumn = 1 + mut->End() + mut->LengthDiff();

    double score;

    bool atBegin = mut->Start() < 3;
    bool atEnd = (mut->End() + 3) > beta_.Columns();

    if (!atBegin && !atEnd) {
        int extendStartCol, extendLength = 2;

        if (mut->Type == MutationType::DELETION) {
            // Future thought: If we revise the semantic of Extra,
            // we can remove the extend and just link alpha and
            // beta directly.
            extendStartCol = mut->Start() - 1;
        } else {
            extendStartCol = mut->Start();
            assert(extendLength <= EXTEND_BUFFER_COLUMNS);
        }
        extendBuffer_.SetDirection(ScaledMatrix::FORWARD);
        recursor_->ExtendAlpha(alpha_, extendStartCol, extendBuffer_, extendLength);
        score = recursor_->LinkAlphaBeta(extendBuffer_, extendLength, beta_, betaLinkCol,
                                         absoluteLinkColumn) +
                alpha_.GetLogProdScales(0, extendStartCol);
    } else if (!atBegin && atEnd) {
        //
        // Extend alpha to end
        //
        size_t extendStartCol = mut->Start() - 1;
        assert(recursor_->tpl_->Length() + 1 > extendStartCol);
        size_t extendLength = recursor_->tpl_->Length() - extendStartCol + 1;

        extendBuffer_.SetDirection(ScaledMatrix::FORWARD);
        recursor_->ExtendAlpha(alpha_, extendStartCol, extendBuffer_, extendLength);
        score = std::log(extendBuffer_(recursor_->read_.Length(), extendLength - 1)) +
                alpha_.GetLogProdScales(0, extendStartCol) +
                extendBuffer_.GetLogProdScales(0, extendLength);
    } else if (atBegin && !atEnd) {
        // If the mutation occurs at positions 0 - 2
        size_t extendLastCol = mut->End();
        // We duplicate this math inside the function
        size_t extendLength = 1 + mut->End() + mut->LengthDiff();

        extendBuffer_.SetDirection(ScaledMatrix::REVERSE);
        recursor_->ExtendBeta(beta_, extendLastCol, extendBuffer_, mut->LengthDiff());
        score = std::log(extendBuffer_(0, 0)) +
                beta_.GetLogProdScales(extendLastCol + 1, beta_.Columns()) +
                extendBuffer_.GetLogProdScales(0, extendLength);
    } else {
        assert(atBegin && atEnd);
        /* This should basically never happen...
           and is a total disaster if it does.  The basic idea is that
           FillAlpha and FillBeta use the real "template" while we test
           mutations using "virtual" template positions and the Extend/Link
           methods.  Trying to call FillAlpha to calculate the likelihood of a
        virtual
           mutation is therefore going to fail, as it calculates using the
           "real" template.
        throw TooSmallTemplateException();
         */

        //
        // Just do the whole fill
        //
        ScaledMatrix alphaP(recursor_->read_.Length() + 1, recursor_->tpl_->Length() + 1,
                            ScaledMatrix::FORWARD);
        recursor_->FillAlpha(ScaledMatrix::Null(), alphaP);
        score = std::log(alphaP(recursor_->read_.Length(), recursor_->tpl_->Length())) +
                alphaP.GetLogProdScales();
    }

    // reset the virtual mutation
    recursor_->tpl_->Reset();

    return score + recursor_->UndoCounterWeights(recursor_->read_.Length());
}

double EvaluatorImpl::LL() const
{
    return std::log(beta_(0, 0)) + beta_.GetLogProdScales() +
           recursor_->UndoCounterWeights(recursor_->read_.Length());
}

std::pair<double, double> EvaluatorImpl::NormalParameters() const
{
    return recursor_->tpl_->NormalParameters();
}

double EvaluatorImpl::ZScore() const
{
    double mean, var;
    std::tie(mean, var) = NormalParameters();
    return (LL() - mean) / std::sqrt(var);
}

inline void EvaluatorImpl::Recalculate()
{
    size_t I = recursor_->read_.Length() + 1;
    size_t J = recursor_->tpl_->Length() + 1;
    alpha_.Reset(I, J);
    beta_.Reset(I, J);
    extendBuffer_.Reset(I, EXTEND_BUFFER_COLUMNS);
    recursor_->FillAlphaBeta(alpha_, beta_);
}

bool EvaluatorImpl::ApplyMutation(const Mutation& mut)
{
    if (recursor_->tpl_->ApplyMutation(mut)) {
        Recalculate();
        return true;
    }
    return false;
}

bool EvaluatorImpl::ApplyMutations(std::vector<Mutation>* muts)
{
    if (recursor_->tpl_->ApplyMutations(muts)) {
        Recalculate();
        return true;
    }
    return false;
}

/*
Matrix* AlphaMatrix();
Matrix* BetaMatrix();
*/

}  // namespace Consensus
}  // namespace PacBio
