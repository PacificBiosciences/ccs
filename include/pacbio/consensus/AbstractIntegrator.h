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

#pragma once

#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <memory>
#include <numeric>

#include <pacbio/consensus/Evaluator.h>
#include <pacbio/consensus/Mutation.h>
#include <pacbio/data/State.h>
#include <pacbio/exception/StateError.h>

namespace PacBio {
namespace Consensus {
// Forward decl
class AbstractMatrix;

/// Contains user-provided filtering information for the Evaluators.
struct IntegratorConfig
{
    double MinZScore;
    double ScoreDiff;

    IntegratorConfig(double minZScore = -3.4, double scoreDiff = 12.5);
};

/// At its core, this class holds a vector of Evaluators and provides helper
/// functions to execute certain actions on each Evaluator.
class AbstractIntegrator
{
public:
    virtual ~AbstractIntegrator();

    virtual size_t TemplateLength() const = 0;

    virtual char operator[](size_t i) const = 0;
    virtual operator std::string() const = 0;

    /// This method throws InvalidEvaluatorException, every time the likelihood
    /// can't be computed for one Evaluator; the respective Evaluator will be invalidated.
    /// You MUST recompute the LLs for all your mutations of interest, as the
    /// number of active Evaluators changed.
    virtual double LL(const Mutation& mut);
    virtual double LL() const;

    double AvgZScore() const;
    std::vector<double> ZScores() const;

    std::vector<std::pair<double, double>> NormalParameters() const;

    virtual void ApplyMutation(const Mutation& mut) = 0;
    virtual void ApplyMutations(std::vector<Mutation>* muts) = 0;

    virtual PacBio::Data::State AddRead(const PacBio::Data::MappedRead& read) = 0;

    /// Given a Mutation of interest, returns a vector of LLs,
    /// one LL per active Evaluator; invalid Evaluators are omitted.
    ///
    /// This method throws InvalidEvaluatorException, every time the likelihood
    /// can't be computed for one Evaluator; the respective Evaluator will be invalidated.
    /// You MUST recompute the LLs for all your mutations of interest, as the
    /// number of active Evaluators changed.
    std::vector<double> LLs(const Mutation& mut);
    /// Return the LL for each Evaluator, even invalid ones.
    /// DO NOT use this in production code, only for debugging purposes.
    std::vector<double> LLs() const;
    /// For each Evaluator, returns the read name.
    std::vector<std::string> ReadNames() const;
    /// Returns the number of flip flop events for each Evaluator.
    std::vector<int> NumFlipFlops() const;
    /// Returns the maximal number of flip flop events of all Evaluators.
    int MaxNumFlipFlops() const;
    /// Computes the ratio of populated cells in the alpha matrix for each
    /// Evaluator and returns the maximal ratio.
    float MaxAlphaPopulated() const;
    /// Computes the ratio of populated cells in the beta matrix for each
    /// Evaluator and returns the maximal ratio.
    float MaxBetaPopulated() const;
    /// Returns the state of each Evaluator.
    std::vector<PacBio::Data::State> States() const;
    /// Returns the strand of each Evaluator.
    std::vector<PacBio::Data::StrandType> StrandTypes() const;

    /// Returns read-only access to Evaluator idx.
    const Evaluator& GetEvaluator(size_t idx) const;

public:
    // Abstract matrix access for SWIG and diagnostics
    const AbstractMatrix& Alpha(size_t idx) const;
    const AbstractMatrix& Beta(size_t idx) const;

protected:
    Mutation ReverseComplement(const Mutation& mut) const;

    AbstractIntegrator(const IntegratorConfig& cfg);

    // move constructor
    AbstractIntegrator(AbstractIntegrator&&);

    PacBio::Data::State AddRead(std::unique_ptr<AbstractTemplate>&& tpl,
                                const PacBio::Data::MappedRead& read);

    IntegratorConfig cfg_;
    std::vector<Evaluator> evals_;

private:
    /// Extract a feature vector from a vector of Evaluators for non-const functions.
    template <typename T>
    inline std::vector<T> TransformEvaluators(std::function<T(Evaluator&)> functor)
    {
        // TODO(atoepfer) How can we use const_cast to convert this?
        std::vector<T> vec;
        vec.reserve(evals_.size());
        std::transform(evals_.begin(), evals_.end(), std::back_inserter(vec), functor);
        return vec;
    }

    /// Extract a feature vector from a vector of Evaluators for const functions.
    template <typename T>
    inline std::vector<T> TransformEvaluators(std::function<T(const Evaluator&)> functor) const
    {
        std::vector<T> vec;
        vec.reserve(evals_.size());
        std::transform(evals_.begin(), evals_.end(), std::back_inserter(vec), functor);
        return vec;
    }
};

/// Helper function to get maximal number from a vector.
template <typename T>
inline T MaxElement(const std::vector<T>& in)
{
    return *std::max_element(in.cbegin(), in.end());
}

}  // namespace Consensus
}  // namespace PacBio
