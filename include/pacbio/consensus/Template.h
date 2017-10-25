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

#include <pacbio/consensus/ModelConfig.h>
#include <pacbio/consensus/Mutation.h>
#include <pacbio/data/Read.h>
#include <pacbio/exception/StateError.h>

namespace PacBio {
namespace Consensus {

// fwd decl
class MutatedTemplate;
class AbstractRecursor;
class ScaledMatrix;

// AbstractTemplate defines the API for representing some provisional
// template or consensus, which need to enable both adding data to
// and updating the underlying sequence
class AbstractTemplate
{
public:
    virtual ~AbstractTemplate();

    size_t Start() const { return start_; }
    virtual size_t Length() const = 0;
    virtual const TemplatePosition& operator[](size_t i) const = 0;

    operator std::string() const;

    // Get a View over the Template (for mutation testing purposes)
    boost::optional<MutatedTemplate> Mutate(const Mutation& m);

    // actually apply mutations
    virtual bool ApplyMutation(const Mutation& mut);
    virtual bool ApplyMutations(std::vector<Mutation>* muts);

    // access model configuration
    virtual std::unique_ptr<AbstractRecursor> CreateRecursor(const PacBio::Data::MappedRead& mr,
                                                             double scoreDiff) const = 0;

    virtual double ExpectedLLForEmission(MoveType move, const NCBI4na prev, const NCBI4na curr,
                                         MomentType moment) const = 0;

    std::pair<double, double> NormalParameters() const;

protected:
    AbstractTemplate(size_t start, size_t end, bool pinStart, bool pinEnd);

    bool InRange(size_t start, size_t end) const;

    virtual const ModelConfig& Config() const = 0;

    size_t start_;
    size_t end_;
    bool pinStart_;
    bool pinEnd_;

private:
    std::pair<double, double> SiteNormalParameters(size_t i) const;

    friend class MutatedTemplate;
};

std::ostream& operator<<(std::ostream&, const AbstractTemplate&);

// Represent a Template sequence as induced by some particular configuration
// of chemistry and model training
class Template : public AbstractTemplate
{
public:
    Template(Template&&) = default;
    Template(const std::string& tpl, std::unique_ptr<ModelConfig>&& cfg);
    Template(const std::string& tpl, std::unique_ptr<ModelConfig>&& cfg, size_t start, size_t end,
             bool pinStart, bool pinEnd);
    size_t Length() const override;
    const TemplatePosition& operator[](size_t i) const override;

    bool ApplyMutation(const Mutation& mut) override;

    std::unique_ptr<AbstractRecursor> CreateRecursor(const PacBio::Data::MappedRead& mr,
                                                     double scoreDiff) const override;

    double ExpectedLLForEmission(MoveType move, const NCBI4na prev, const NCBI4na curr,
                                 MomentType moment) const override;

protected:
    const ModelConfig& Config() const override { return *cfg_; }

private:
    std::unique_ptr<ModelConfig> cfg_;
    std::vector<TemplatePosition> tpl_;
};

// A View projected from some template, allowing for the analysis of a
// hypothetical mutation without modifying the underlying Template,
// which can now be kept const
class MutatedTemplate : public AbstractTemplate
{
public:
    MutatedTemplate(const AbstractTemplate& master, const Mutation& m);

    MutatedTemplate(const MutatedTemplate&) = default;
    MutatedTemplate(MutatedTemplate&&) = default;

    size_t Length() const override;
    MutationType Type() const;
    size_t MutationStart() const;
    size_t MutationEnd() const;
    int LengthDiff() const;
    const TemplatePosition& operator[](size_t i) const override;

    bool ApplyMutation(const Mutation& mut) override;

    std::unique_ptr<AbstractRecursor> CreateRecursor(const PacBio::Data::MappedRead& mr,
                                                     double scoreDiff) const override;

    double ExpectedLLForEmission(MoveType move, const NCBI4na prev, const NCBI4na curr,
                                 MomentType moment) const override;

protected:
    const ModelConfig& Config() const override { return master_.Config(); }

private:
    const AbstractTemplate& master_;
    const Mutation mut_;
    const size_t mutStart_;
    const int mutOff_;
    std::vector<TemplatePosition> mutTpl_;  // params for the context start at the mutation
};

// this needs to be here because the unique_ptr deleter for AbstractRecursor must know its size
class AbstractRecursor
{
protected:
    typedef ScaledMatrix M;

public:
    AbstractRecursor(PacBio::Data::MappedRead mr, double scoreDiff);
    virtual ~AbstractRecursor() {}
    virtual size_t FillAlphaBeta(const AbstractTemplate& tpl, M& alpha, M& beta,
                                 double tol) const = 0;
    virtual void FillAlpha(const AbstractTemplate& tpl, const M& guide, M& alpha) const = 0;
    virtual void FillBeta(const AbstractTemplate& tpl, const M& guide, M& beta) const = 0;
    virtual double LinkAlphaBeta(const AbstractTemplate& tpl, const M& alpha, size_t alphaColumn,
                                 const M& beta, size_t betaColumn, size_t absoluteColumn) const = 0;
    virtual void ExtendAlpha(const AbstractTemplate& tpl, const M& alpha, size_t beginColumn,
                             M& ext, size_t numExtColumns = 2) const = 0;
    virtual void ExtendBeta(const AbstractTemplate& tpl, const M& beta, size_t endColumn, M& ext,
                            int lengthDiff = 0) const = 0;
    virtual double UndoCounterWeights(size_t nEmissions) const = 0;

public:
    PacBio::Data::MappedRead read_;
    const double scoreDiff_;  // reciprocal of "natural scale"
};

}  // namespace Consensus
}  // namespace PacBio
