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
#include <set>

#include <pacbio/consensus/Evaluator.h>
#include <pacbio/consensus/Exceptions.h>
#include <pacbio/consensus/Mutation.h>
#include <pacbio/consensus/State.h>

namespace PacBio {
namespace Consensus {

struct IntegratorConfig
{
    double MinZScore;
    double ScoreDiff;

    IntegratorConfig(double minZScore = -3.5, double scoreDiff = 12.5);
};

inline std::ostream& operator<<(std::ostream& os, State result)
{
    os << StateName[static_cast<size_t>(result)];
    return os;
}

std::set<std::string> SupportedChemistries();

class AbstractIntegrator
{
public:
    virtual ~AbstractIntegrator();

    virtual size_t TemplateLength() const = 0;

    virtual char operator[](size_t i) const = 0;
    virtual operator std::string() const = 0;

    virtual double LL(const Mutation& mut);
    virtual double LL() const;

    double AvgZScore() const;
    std::vector<double> ZScores() const;

    std::vector<std::pair<double, double>> NormalParameters() const;

    virtual void ApplyMutation(const Mutation& mut) = 0;
    virtual void ApplyMutations(std::vector<Mutation>* muts) = 0;

    virtual State AddRead(const MappedRead& read) = 0;

    // For debugging purposes
    // (Note that these include results include all evaluators, even the inactive ones)
    std::vector<double> LLs(const Mutation& mut);
    std::vector<double> LLs() const;
    std::vector<std::string> ReadNames() const;

    std::vector<int> NumFlipFlops() const;
    int MaxNumFlipFlops() const;

    std::vector<float> AlphaPopulated() const;
    float MaxAlphaPopulated() const;

    std::vector<float> BetaPopulated() const;
    float MaxBetaPopulated() const;

    std::vector<State> States() const;
    std::vector<StrandType> StrandTypes() const;

    // TODO(atoepfer) Does anyone have a clue if we can make one function out
    //                of those two?
    template <typename T>
    inline std::vector<T> TransformEvaluators(std::function<T(Evaluator&)> functor)
    {
        std::vector<T> vec;
        vec.reserve(evals_.size());
        std::transform(evals_.begin(), evals_.end(), std::back_inserter(vec), functor);
        return vec;
    }

    template <typename T>
    inline std::vector<T> TransformEvaluators(std::function<T(const Evaluator&)> functor) const
    {
        std::vector<T> vec;
        vec.reserve(evals_.size());
        std::transform(evals_.begin(), evals_.end(), std::back_inserter(vec), functor);
        return vec;
    }

    template <typename T>
    inline T MaxElement(const std::vector<T>& in) const
    {
        return *std::max_element(in.cbegin(), in.end());
    }

protected:
    Mutation ReverseComplement(const Mutation& mut) const;

    AbstractIntegrator(const IntegratorConfig& cfg);

    // move constructor
    AbstractIntegrator(AbstractIntegrator&&);

    State AddRead(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& read);

    IntegratorConfig cfg_;
    std::vector<Evaluator> evals_;

private:
    inline double AccumulateNoInf(std::vector<double> input) const
    {
        const auto AddNoInf = [](double a, double b) { return a + (std::isinf(b) ? 0.0 : b); };
        return std::accumulate(input.cbegin(), input.cend(), 0.0, AddNoInf);
    }
};

}  // namespace Consensus
}  // namespace PacBio
