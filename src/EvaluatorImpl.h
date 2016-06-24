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

#include <pacbio/consensus/Evaluator.h>
#include <pacbio/consensus/Read.h>
#include <pacbio/consensus/Template.h>

#include "Recursor.h"
#include "matrix/ScaledMatrix.h"

namespace PacBio {
namespace Consensus {

class EvaluatorImpl
{
public:
    EvaluatorImpl(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr,
                  double scoreDiff = 12.5);

    std::string ReadName() const;

    double LL(const Mutation& mut);
    double LL() const;

    // TODO: Comments are nice!  Explain what this is about---ZScore calculation?
    std::pair<double, double> NormalParameters() const;

    double ZScore() const;

    bool ApplyMutation(const Mutation& mut);
    bool ApplyMutations(std::vector<Mutation>* muts);

private:
    void Recalculate();

private:
    std::unique_ptr<AbstractRecursor>
        recursor_;  // TODO: does this need to be a pointer?  is it always non-null?
                    // are we making it a UP just so we can do a fwd decl and still have
                    // RAII semantics?
    ScaledMatrix alpha_;
    ScaledMatrix beta_;
    ScaledMatrix extendBuffer_;

    friend class Evaluator;
};

}  // namespace Consensus
}  // namespace PacBio
