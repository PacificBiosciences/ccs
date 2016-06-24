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

#include <pacbio/consensus/Read.h>
#include <pacbio/consensus/Template.h>

namespace PacBio {
namespace Consensus {

//
// TODO: explain the legal state transitions.  What does it mean to be in each state?
// Can a read "come back" form being DISABLED to being VALID (I hope not!) or from one
// of the other states?
// TODO: Get rid of NULL_TEMPLATE state---this should be handled by ReadScoresMutation logic
enum struct EvaluatorState : uint8_t
{
    VALID,
    ALPHA_BETA_MISMATCH,
    POOR_ZSCORE,
    NULL_TEMPLATE,
    DISABLED
};

// forward declaration
class EvaluatorImpl;

class Evaluator
{
public:
    Evaluator(EvaluatorState);
    Evaluator(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr, double minZScore,
              double scoreDiff);

    // move constructors
    Evaluator(Evaluator&&);
    Evaluator& operator=(Evaluator&&);

    ~Evaluator();

    size_t Length() const;  // TODO: is this used anywhere?  If not, delete it.
    StrandEnum Strand() const;

    operator bool() const;
    operator std::string() const;
    std::string ReadName() const;

    double LL(const Mutation& mut);
    double LL() const;

    std::pair<double, double> NormalParameters() const;

    double ZScore() const;

    bool ApplyMutation(const Mutation& mut);
    bool ApplyMutations(std::vector<Mutation>* muts);

    EvaluatorState Status() const;
    void Release();

private:
    void CheckInvariants();

private:
    std::unique_ptr<EvaluatorImpl> impl_;
    EvaluatorState state_;
};

}  // namespace Consensus
}  // namespace PacBio
