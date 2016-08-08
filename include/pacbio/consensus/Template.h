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
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <boost/optional.hpp>

#include <pacbio/exception/StateError.h>
#include <pacbio/consensus/ModelConfig.h>
#include <pacbio/consensus/Mutation.h>
#include <pacbio/data/Read.h>

namespace PacBio {
namespace Consensus {

// fwd decl
class AbstractRecursor;
class ScaledMatrix;

class AbstractTemplate
{
public:
    virtual ~AbstractTemplate();

    virtual size_t Length() const = 0;
    virtual const TemplatePosition& operator[](size_t i) const = 0;

    operator std::string() const;

    // virtual mutations (for mutation testing purposes)
    virtual bool IsMutated() const = 0;
    virtual boost::optional<Mutation> Mutate(const Mutation& m) = 0;
    virtual void Reset() = 0;

    // actually apply mutations
    virtual bool ApplyMutation(const Mutation& mut);
    virtual bool ApplyMutations(std::vector<Mutation>* muts);

    // access model configuration
    virtual std::unique_ptr<AbstractRecursor> CreateRecursor(
        std::unique_ptr<AbstractTemplate>&& tpl, const PacBio::Data::MappedRead& mr, double scoreDiff) const = 0;
    virtual double ExpectedLLForEmission(MoveType move, uint8_t prev, uint8_t curr,
                                         MomentType moment) const = 0;

    std::pair<double, double> NormalParameters() const;

    // a sad but necessary release valve for MonoMolecularIntegrator Length()
    size_t TrueLength() const;

protected:
    AbstractTemplate(size_t start, size_t end, bool pinStart, bool pinEnd);

    bool InRange(size_t start, size_t end) const;

    size_t start_;
    size_t end_;
    bool pinStart_;
    bool pinEnd_;

private:
    std::pair<double, double> SiteNormalParameters(size_t i) const;
};

std::ostream& operator<<(std::ostream&, const AbstractTemplate&);

class Template : public AbstractTemplate
{
public:
    Template(const std::string& tpl, std::unique_ptr<ModelConfig>&& cfg);
    Template(const std::string& tpl, std::unique_ptr<ModelConfig>&& cfg, size_t start, size_t end,
             bool pinStart, bool pinEnd);

    size_t Length() const override;
    const TemplatePosition& operator[](size_t i) const override;

    bool IsMutated() const override;
    boost::optional<Mutation> Mutate(const Mutation& mut) override;
    void Reset() override;

    bool ApplyMutation(const Mutation& mut) override;

    std::unique_ptr<AbstractRecursor> CreateRecursor(std::unique_ptr<AbstractTemplate>&& tpl,
                                                     const PacBio::Data::MappedRead& mr,
                                                     double scoreDiff) const override;

    double ExpectedLLForEmission(MoveType move, uint8_t prev, uint8_t curr,
                                 MomentType moment) const override;

private:
    std::unique_ptr<ModelConfig> cfg_;
    std::vector<TemplatePosition> tpl_;
    bool mutated_;
    size_t mutStart_;
    size_t mutEnd_;
    int mutOff_;
    TemplatePosition mutTpl_[2];

    friend class VirtualTemplate;
};

class VirtualTemplate : public AbstractTemplate
{
public:
    VirtualTemplate(const Template& master, size_t start, size_t end, bool pinStart, bool pinEnd);

    size_t Length() const override;
    const TemplatePosition& operator[](size_t i) const override;

    bool IsMutated() const override;
    boost::optional<Mutation> Mutate(const Mutation&) override;
    void Reset() override {}
    bool ApplyMutation(const Mutation& mut) override;

    std::unique_ptr<AbstractRecursor> CreateRecursor(std::unique_ptr<AbstractTemplate>&& tpl,
                                                     const PacBio::Data::MappedRead& mr,
                                                     double scoreDiff) const override;

    double ExpectedLLForEmission(MoveType move, uint8_t prev, uint8_t curr,
                                 MomentType moment) const override;

private:
    Template const& master_;
};

// this needs to be here because the unique_ptr deleter for AbstractRecursor must know its size
class AbstractRecursor
{
protected:
    typedef ScaledMatrix M;

public:
    AbstractRecursor(std::unique_ptr<AbstractTemplate>&& tpl, const PacBio::Data::MappedRead& mr,
                     double scoreDiff);
    virtual ~AbstractRecursor() {}
    virtual size_t FillAlphaBeta(M& alpha, M& beta) const = 0;
    virtual void FillAlpha(const M& guide, M& alpha) const = 0;
    virtual void FillBeta(const M& guide, M& beta) const = 0;
    virtual double LinkAlphaBeta(const M& alpha, size_t alphaColumn, const M& beta,
                                 size_t betaColumn, size_t absoluteColumn) const = 0;
    virtual void ExtendAlpha(const M& alpha, size_t beginColumn, M& ext,
                             size_t numExtColumns = 2) const = 0;
    virtual void ExtendBeta(const M& beta, size_t endColumn, M& ext, int lengthDiff = 0) const = 0;
    virtual double UndoCounterWeights(size_t nEmissions) const = 0;

public:
    std::unique_ptr<AbstractTemplate> tpl_;
    PacBio::Data::MappedRead read_;

protected:
    double scoreDiff_;  // reciprocal of "natural scale"
};

}  // namespace Consensus
}  // namespace PacBio
