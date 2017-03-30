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

#include <cassert>
#include <cmath>
#include <limits>
#include <utility>

#include <pacbio/consensus/Integrator.h>
#include <pacbio/data/Sequence.h>

#include "ModelFactory.h"

using namespace PacBio::Data;
using namespace PacBio::Exception;

namespace PacBio {
namespace Consensus {

Integrator::Integrator(const std::string& tpl, const IntegratorConfig& cfg)
    : AbstractIntegrator(cfg), fwdTpl_{tpl}, revTpl_{::PacBio::Data::ReverseComplement(tpl)}
{
}

State Integrator::AddRead(const PacBio::Data::MappedRead& read)
{
    try {
        return AbstractIntegrator::AddRead(GetTemplate(read), read);
    } catch (const TemplateTooSmall& e) {
        return State::TEMPLATE_TOO_SMALL;
    }
}

size_t Integrator::TemplateLength() const { return fwdTpl_.length(); }

char Integrator::operator[](const size_t i) const { return fwdTpl_[i]; }

Integrator::operator std::string() const { return fwdTpl_; }

void Integrator::ApplyMutation(const Mutation& fwdMut)
{
    const Mutation revMut(ReverseComplement(fwdMut));

    std::vector<Mutation> fwdMuts = {fwdMut};
    std::vector<Mutation> revMuts = {revMut};

    fwdTpl_ = ::PacBio::Consensus::ApplyMutations(fwdTpl_, &fwdMuts);
    revTpl_ = ::PacBio::Consensus::ApplyMutations(revTpl_, &revMuts);

    for (auto& eval : evals_) {
        if (eval.Strand() == StrandType::FORWARD)
            eval.ApplyMutation(fwdMut);
        else if (eval.Strand() == StrandType::REVERSE)
            eval.ApplyMutation(revMut);
    }

    assert(fwdTpl_.length() == revTpl_.length());
    assert(fwdTpl_ == ::PacBio::Data::ReverseComplement(revTpl_));
}

void Integrator::ApplyMutations(std::vector<Mutation>* fwdMuts)
{
    std::vector<Mutation> revMuts;

    for (auto it = fwdMuts->crbegin(); it != fwdMuts->crend(); ++it)
        revMuts.emplace_back(ReverseComplement(*it));

    fwdTpl_ = ::PacBio::Consensus::ApplyMutations(fwdTpl_, fwdMuts);
    revTpl_ = ::PacBio::Consensus::ApplyMutations(revTpl_, &revMuts);

    for (auto& eval : evals_) {
        if (eval.Strand() == StrandType::FORWARD)
            eval.ApplyMutations(fwdMuts);
        else if (eval.Strand() == StrandType::REVERSE)
            eval.ApplyMutations(&revMuts);
    }

    assert(fwdTpl_.length() == revTpl_.length());
    assert(fwdTpl_ == ::PacBio::Data::ReverseComplement(revTpl_));
}

std::unique_ptr<AbstractTemplate> Integrator::GetTemplate(const PacBio::Data::MappedRead& read)
{
    const size_t len = read.TemplateEnd - read.TemplateStart;

    if (read.Strand == StrandType::FORWARD) {
        const size_t start = read.TemplateStart;
        const size_t end = read.TemplateEnd;

        return std::unique_ptr<AbstractTemplate>(new Template(fwdTpl_.substr(start, len),
                                                              ModelFactory::Create(read), start,
                                                              end, read.PinStart, read.PinEnd));
    } else if (read.Strand == StrandType::REVERSE) {
        const size_t start = revTpl_.size() - read.TemplateEnd;
        const size_t end = revTpl_.size() - read.TemplateStart;

        return std::unique_ptr<AbstractTemplate>(new Template(revTpl_.substr(start, len),
                                                              ModelFactory::Create(read), start,
                                                              end, read.PinEnd, read.PinStart));
    }

    throw std::invalid_argument("read is unmapped!");
}

}  // namespace Consensus
}  // namespace PacBio
