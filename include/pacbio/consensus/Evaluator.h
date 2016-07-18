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

#include <memory>
#include <utility>
#include <vector>

#include <pacbio/data/Read.h>
#include <pacbio/data/State.h>
#include <pacbio/consensus/Template.h>

namespace PacBio {
namespace Consensus {
    
// forward declaration
class EvaluatorImpl;

class Evaluator
{
public:
    Evaluator() = delete;
    Evaluator(PacBio::Data::State);
    Evaluator(std::unique_ptr<AbstractTemplate>&& tpl, const PacBio::Data::MappedRead& mr, double minZScore,
              double scoreDiff);

    // copying is verboten
    Evaluator(const Evaluator&) = delete;
    Evaluator& operator=(const Evaluator&) = delete;

    // move constructor
    Evaluator(Evaluator&&);
    // move assign operator
    Evaluator& operator=(Evaluator&&);

    ~Evaluator();

    size_t Length() const;  // TODO: is this used anywhere?  If not, delete it.
    PacBio::Data::StrandType Strand() const;

    operator bool() const { return IsValid(); }
    operator std::string() const;
    std::string ReadName() const;

    double LL(const Mutation& mut);
    double LL() const;

    std::pair<double, double> NormalParameters() const;

    double ZScore() const;

    bool ApplyMutation(const Mutation& mut);
    bool ApplyMutations(std::vector<Mutation>* muts);

    PacBio::Data::State Status() const { return curState_; }
    int NumFlipFlops() const;
    float AlphaPopulated() const;
    float BetaPopulated() const;

    void Release();

private:
    void CheckZScore(const double minZScore, const std::string& model);

    bool IsValid() const { return curState_ == PacBio::Data::State::VALID; }
    void Status(PacBio::Data::State nextState);

private:
    std::unique_ptr<EvaluatorImpl> impl_;
    PacBio::Data::State curState_;
};

}  // namespace Consensus
}  // namespace PacBio
