// Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.
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

// Author: Patrick Marks and David Alexander

#include <ConsensusCore/Quiver/MutationScorer.hpp>

#include <ConsensusCore/Edna/EdnaEvaluator.hpp>
#include <ConsensusCore/Matrix/DenseMatrix.hpp>
#include <ConsensusCore/Matrix/SparseMatrix.hpp>
#include <ConsensusCore/Quiver/QvEvaluator.hpp>
#include <ConsensusCore/Quiver/SimpleRecursor.hpp>
#include <ConsensusCore/Quiver/SseRecursor.hpp>
#include <ConsensusCore/Mutation.hpp>

#include <string>

#define EXTEND_BUFFER_COLUMNS 8

namespace ConsensusCore
{
    template<typename R>
    MutationScorer<R>::MutationScorer(const EvaluatorType& evaluator, const R& recursor)
        throw(AlphaBetaMismatchException)
        : evaluator_(new EvaluatorType(evaluator)),
          recursor_(new R(recursor))
    {
        try {
            // Allocate alpha and beta
            alpha_ = new MatrixType(evaluator.ReadLength() + 1,
                                    evaluator.TemplateLength() + 1);
            beta_ = new MatrixType(evaluator.ReadLength() + 1,
                                   evaluator.TemplateLength() + 1);
            // Buffer where we extend into
            extendBuffer_ = new MatrixType(evaluator.ReadLength() + 1, EXTEND_BUFFER_COLUMNS);
            // Initial alpha and beta
            numFlipFlops_ = recursor.FillAlphaBeta(*evaluator_, *alpha_, *beta_);
        }
        catch(AlphaBetaMismatchException e) {
            delete alpha_;
            delete beta_;
            delete extendBuffer_;
            delete recursor_;
            throw;
        }
    }

    template<typename R>
    MutationScorer<R>::MutationScorer(const MutationScorer<R>& other)
    {
        evaluator_ = new EvaluatorType(*other.evaluator_);
        recursor_ = new R(*other.recursor_);

        // Copy alpha and beta
        alpha_ = new MatrixType(*other.alpha_);
        beta_ = new MatrixType(*other.beta_);
        // Buffer where we extend into
        extendBuffer_ = new MatrixType(*other.extendBuffer_);
        numFlipFlops_ = other.numFlipFlops_;
    }

    template<typename R>
    float
    MutationScorer<R>::Score() const
    {
        return (*beta_)(0, 0);
    }

    template<typename R> std::string
    MutationScorer<R>::Template() const
    {
        return evaluator_->Template();
    }

    template<typename R>
    void MutationScorer<R>::Template(std::string tpl)
        throw(AlphaBetaMismatchException)
    {
        delete alpha_;
        delete beta_;
        evaluator_->Template(tpl);
        alpha_ = new MatrixType(evaluator_->ReadLength() + 1,
                                evaluator_->TemplateLength() + 1);
        beta_  = new MatrixType(evaluator_->ReadLength() + 1,
                                evaluator_->TemplateLength() + 1);
        recursor_->FillAlphaBeta(*evaluator_, *alpha_, *beta_);
    }

    template<typename R>
    const typename R::MatrixType* MutationScorer<R>::Alpha() const
    {
        return alpha_;
    }

    template<typename R>
    const typename R::MatrixType* MutationScorer<R>::Beta() const
    {
        return beta_;
    }

    template<typename R>
    const typename R::EvaluatorType* MutationScorer<R>::Evaluator() const
    {
        return evaluator_;
    }

    template<typename R>
    const PairwiseAlignment* MutationScorer<R>::Alignment() const
    {
        return recursor_->Alignment(*evaluator_, *alpha_);
    }

    template<typename R>
    float
    MutationScorer<R>::ScoreMutation(const Mutation& m) const
    {
        int betaLinkCol = 1 + m.End();
        int absoluteLinkColumn = 1 + m.End() + m.LengthDiff();
        std::string oldTpl = evaluator_->Template();
        std::string newTpl = ApplyMutation(m, oldTpl);
        float score;

        bool atBegin = (m.Start() < 3);
        bool atEnd   = (m.End() > (int)oldTpl.length() - 2);

        if (!atBegin && !atEnd)
        {
            // Install mutated template
            evaluator_->Template(newTpl);

            int extendStartCol, extendLength;

            if (m.Type() == DELETION)
            {
                // Future thought: If we revise the semantic of Extra,
                // we can remove the extend and just link alpha and
                // beta directly.
                extendStartCol = m.Start() - 1;
                extendLength = 2;
            }
            else
            {
                extendStartCol = m.Start();
                extendLength   = 1 + m.NewBases().length();
                assert(extendLength <= EXTEND_BUFFER_COLUMNS);
            }

            recursor_->ExtendAlpha(*evaluator_, *alpha_,
                                   extendStartCol, *extendBuffer_, extendLength);
            score = recursor_->LinkAlphaBeta(*evaluator_,
                                             *extendBuffer_, extendLength,
                                             *beta_, betaLinkCol,
                                             absoluteLinkColumn);
        }
        else if (!atBegin && atEnd)
        {
            //
            // Extend alpha to end
            //
            evaluator_->Template(newTpl);

            int extendStartCol = m.Start() - 1;
            int extendLength = newTpl.length() - extendStartCol + 1;

            recursor_->ExtendAlpha(*evaluator_, *alpha_,
                                   extendStartCol, *extendBuffer_, extendLength);
            score = (*extendBuffer_)(evaluator_->ReadLength(), extendLength - 1);

            // if (fabs(score - Score()) > 50) {
            //     // FIXME!  This happens on fluidigm amplicons, figure out why
            //     Breakpoint();
            // }
        }
        else if (atBegin && !atEnd)
        {
            //
            // Extend beta back
            //
            evaluator_->Template(newTpl);

            int extendLastCol = m.End();
            int extendLength = m.End() + m.LengthDiff() + 1;

            recursor_->ExtendBeta(*evaluator_, *beta_,
                                  extendLastCol, *extendBuffer_, extendLength,
                                  m.LengthDiff());
            score = (*extendBuffer_)(0, 0);
        }
        else
        {
            assert(atBegin && atEnd);
            //
            // Just do the whole fill
            //
            MatrixType alphaP(evaluator_->ReadLength() + 1,
                              newTpl.length() + 1);
            evaluator_->Template(newTpl);
            recursor_->FillAlpha(*evaluator_, MatrixType::Null(), alphaP);
            score = alphaP(evaluator_->ReadLength(), newTpl.length());
        }

        // Restore the original template.
        evaluator_->Template(oldTpl);

        // if (fabs(score - Score()) > 50) { Breakpoint(); }

        return score;
    }


    template<typename R>
    MutationScorer<R>::~MutationScorer()
    {
        delete extendBuffer_;
        delete beta_;
        delete alpha_;
        delete recursor_;
        delete evaluator_;
    }

    template class MutationScorer<SimpleQvRecursor>;
    template class MutationScorer<SseQvRecursor>;
    template class MutationScorer<SparseSimpleQvRecursor>;
    template class MutationScorer<SparseSimpleQvSumProductRecursor>;
    template class MutationScorer<SparseSseQvRecursor>;
    template class MutationScorer<SparseSseQvSumProductRecursor>;
    template class MutationScorer<SparseSseEdnaRecursor>;
}
