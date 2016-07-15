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

#include <cstdint>
#include <functional>
#include <iostream>
#include <memory>
#include <set>

#include <pacbio/consensus/AbstractIntegrator.h>
#include <pacbio/consensus/Evaluator.h>
#include <pacbio/consensus/Exceptions.h>
#include <pacbio/consensus/Mutation.h>
#include <pacbio/consensus/State.h>

namespace PacBio {
namespace Consensus {

class MonoMolecularIntegrator : public AbstractIntegrator
{
public:
    MonoMolecularIntegrator(const std::string& tpl, const IntegratorConfig& cfg, const SNR& snr,
                            const std::string& model);

    // move constructor
    MonoMolecularIntegrator(MonoMolecularIntegrator&&);

    size_t TemplateLength() const;

    char operator[](size_t i) const;
    operator std::string() const;

    double LL(const Mutation& mut);
    inline double LL() const { return AbstractIntegrator::LL(); }
    void ApplyMutation(const Mutation& mut);
    void ApplyMutations(std::vector<Mutation>* muts);

    State AddRead(const MappedRead& read);

protected:
    std::string mdl_;
    SNR snr_;
    Template fwdTpl_;
    Template revTpl_;
};

}  // namespace Consensus
}  // namespace PacBio
