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
#include <pacbio/consensus/Mutation.h>
#include <pacbio/data/Read.h>
#include <pacbio/data/State.h>
#include <pacbio/exception/StateError.h>

namespace PacBio {
namespace Consensus {

/// The MULTI-molecular integrator holds those Evaluators, whose MappedReads
/// belong to the same genomic region, but do not share the same template.
class MultiMolecularIntegrator : public AbstractIntegrator
{
public:
    /// \brief Initialize the MultiMolecularIntegrator.
    ///
    /// \param tpl    The draft template as a string
    /// \param cfg    The configuration used to initialize the AbstractIntegrator.
    MultiMolecularIntegrator(const std::string& tpl, const IntegratorConfig& cfg);

    size_t TemplateLength() const override;

    /// Returns base i of the template
    char operator[](size_t i) const override;
    operator std::string() const override;

    /// Applies a mutation to the template of each Evaluator.
    void ApplyMutation(const Mutation& mut) override;
    /// Applies a vector of murations to the template of each Evaluator.
    void ApplyMutations(std::vector<Mutation>* muts) override;
    /// Encapsulate the read in an Evaluator and stores it.
    PacBio::Data::State AddRead(const PacBio::Data::MappedRead& read) override;

protected:
    std::unique_ptr<AbstractTemplate> GetTemplate(const PacBio::Data::MappedRead& read);

    std::string fwdTpl_;
    std::string revTpl_;

private:
    friend struct std::hash<MultiMolecularIntegrator>;
};

}  // namespace Consensus
}  // namespace PacBio
