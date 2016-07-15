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

#include <pacbio/consensus/MonoMolecularIntegrator.h>
#include <pacbio/consensus/Sequence.h>

#include "ModelFactory.h"

namespace PacBio {
namespace Consensus {

MonoMolecularIntegrator::MonoMolecularIntegrator(const std::string& tpl,
                                                 const IntegratorConfig& cfg, const SNR& snr,
                                                 const std::string& model)
    : AbstractIntegrator(cfg)
    , mdl_{model}
    , snr_{snr}
    , fwdTpl_(tpl, ModelFactory::Create(mdl_, snr_))
    , revTpl_(::PacBio::Consensus::ReverseComplement(tpl), ModelFactory::Create(mdl_, snr_))
{
}

MonoMolecularIntegrator::MonoMolecularIntegrator(MonoMolecularIntegrator&& mmi)
    : AbstractIntegrator(std::move(mmi))
    , mdl_{mmi.mdl_}
    , snr_{mmi.snr_}
    , fwdTpl_{std::move(mmi.fwdTpl_)}
    , revTpl_{std::move(mmi.revTpl_)}
{
}

State MonoMolecularIntegrator::AddRead(const MappedRead& read)
{
    if (read.Model != mdl_) throw std::invalid_argument("invalid model for integrator!");
    if (read.SignalToNoise != snr_) throw std::invalid_argument("invalid SNR for integrator!");

    try {
        if (read.Strand == StrandType::FORWARD)
            return AbstractIntegrator::AddRead(
                std::unique_ptr<AbstractTemplate>(new VirtualTemplate(
                    fwdTpl_, read.TemplateStart, read.TemplateEnd, read.PinStart, read.PinEnd)),
                read);

        else if (read.Strand == StrandType::REVERSE)
            return AbstractIntegrator::AddRead(
                std::unique_ptr<AbstractTemplate>(new VirtualTemplate(
                    revTpl_, TemplateLength() - read.TemplateEnd,
                    TemplateLength() - read.TemplateStart, read.PinEnd, read.PinStart)),
                read);
    } catch (const TemplateTooSmall& e) {
        return State::TEMPLATE_TOO_SMALL;
    }

    throw std::invalid_argument("read is unmapped!");
}

size_t MonoMolecularIntegrator::TemplateLength() const { return fwdTpl_.TrueLength(); }
char MonoMolecularIntegrator::operator[](const size_t i) const { return fwdTpl_[i].Base; }
MonoMolecularIntegrator::operator std::string() const
{
    std::string result;

    result.resize(fwdTpl_.Length());

    for (size_t i = 0; i < fwdTpl_.Length(); ++i)
        result[i] = fwdTpl_[i].Base;

    return result;
}

double MonoMolecularIntegrator::LL(const Mutation& fwdMut)
{
    const Mutation revMut(ReverseComplement(fwdMut));
    fwdTpl_.Mutate(fwdMut);
    revTpl_.Mutate(revMut);
    const double ll = AbstractIntegrator::LL(fwdMut);
    fwdTpl_.Reset();
    revTpl_.Reset();
    return ll;
}

void MonoMolecularIntegrator::ApplyMutation(const Mutation& fwdMut)
{
    const Mutation revMut(ReverseComplement(fwdMut));

    fwdTpl_.ApplyMutation(fwdMut);
    revTpl_.ApplyMutation(revMut);

    for (auto& eval : evals_) {
        if (eval.Strand() == StrandType::FORWARD)
            eval.ApplyMutation(fwdMut);
        else if (eval.Strand() == StrandType::REVERSE)
            eval.ApplyMutation(revMut);
    }

    assert(fwdTpl_.Length() == revTpl_.Length());

#ifndef NDEBUG
    std::string fwd;
    std::string rev;

    for (size_t i = 0; i < TemplateLength(); ++i) {
        fwd.push_back(fwdTpl_[i].Base);
        rev.push_back(revTpl_[i].Base);
    }

#endif

    assert(fwd == ::PacBio::Consensus::ReverseComplement(rev));
}

void MonoMolecularIntegrator::ApplyMutations(std::vector<Mutation>* fwdMuts)
{
    std::vector<Mutation> revMuts;

    for (auto it = fwdMuts->crbegin(); it != fwdMuts->crend(); ++it)
        revMuts.emplace_back(ReverseComplement(*it));

    fwdTpl_.ApplyMutations(fwdMuts);
    revTpl_.ApplyMutations(&revMuts);

    for (auto& eval : evals_) {
        if (eval.Strand() == StrandType::FORWARD)
            eval.ApplyMutations(fwdMuts);
        else if (eval.Strand() == StrandType::REVERSE)
            eval.ApplyMutations(&revMuts);
    }

    assert(fwdTpl_.Length() == revTpl_.Length());

#ifndef NDEBUG
    std::string fwd;
    std::string rev;

    for (size_t i = 0; i < TemplateLength(); ++i) {
        fwd.push_back(fwdTpl_[i].Base);
        rev.push_back(revTpl_[i].Base);
    }
#endif

    assert(fwd == ::PacBio::Consensus::ReverseComplement(rev));
}

}  // namespace Consensus
}  // namespace PacBio
