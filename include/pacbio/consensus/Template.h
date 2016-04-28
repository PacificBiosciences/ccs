
#pragma once

#include <cstdint>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <boost/optional.hpp>

#include <pacbio/consensus/Exceptions.h>
#include <pacbio/consensus/ModelConfig.h>
#include <pacbio/consensus/Mutation.h>
#include <pacbio/consensus/Read.h>

namespace PacBio {
namespace Consensus {

// fwd decl
class AbstractRecursor;
class ScaledMatrix;
struct MappedRead;

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
        std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr, double scoreDiff) const = 0;
    virtual double SubstitutionRate(uint8_t prev, uint8_t curr) const = 0;

    std::pair<double, double> NormalParameters() const;

    // a sad but necessary release valve for MonoMolecularIntegrator Length()
    size_t TrueLength() const;

protected:
    AbstractTemplate(size_t start, size_t end, bool pinStart, bool pinEnd);

    inline bool InRange(size_t start, size_t end) const;

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

    inline size_t Length() const;
    inline const TemplatePosition& operator[](size_t i) const;

    inline bool IsMutated() const;
    boost::optional<Mutation> Mutate(const Mutation& mut);
    void Reset();

    bool ApplyMutation(const Mutation& mut);

    inline std::unique_ptr<AbstractRecursor> CreateRecursor(std::unique_ptr<AbstractTemplate>&& tpl,
                                                            const MappedRead& mr,
                                                            double scoreDiff) const;
    inline double SubstitutionRate(uint8_t prev, uint8_t curr) const;

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

    inline size_t Length() const;
    inline const TemplatePosition& operator[](size_t i) const;

    inline bool IsMutated() const;
    inline boost::optional<Mutation> Mutate(const Mutation&);
    inline void Reset() {}
    bool ApplyMutation(const Mutation& mut);

    inline std::unique_ptr<AbstractRecursor> CreateRecursor(std::unique_ptr<AbstractTemplate>&& tpl,
                                                            const MappedRead& mr,
                                                            double scoreDiff) const;
    inline double SubstitutionRate(uint8_t prev, uint8_t curr) const;

private:
    Template const& master_;
};

// this needs to be here because the unique_ptr deleter for AbstractRecursor must know its size
class AbstractRecursor
{
protected:
    typedef ScaledMatrix M;

public:
    AbstractRecursor(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr,
                     double scoreDiff);
    virtual ~AbstractRecursor() {}
    virtual size_t FillAlphaBeta(M& alpha, M& beta) const throw(AlphaBetaMismatch) = 0;
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
    MappedRead read_;

protected:
    double scoreDiff_;  // reciprocal of "natural scale"
};

// inline function impls
bool AbstractTemplate::InRange(const size_t start, const size_t end) const
{
    if ((pinStart_ || start_ < end) && (pinEnd_ || start < end_)) return true;
    return false;
}

size_t Template::Length() const { return tpl_.size() + mutOff_; }

const TemplatePosition& Template::operator[](size_t i) const
{
    // if no mutation, or everything up to the base before mutStart_, just return
    // what we have
    if (!IsMutated() || i + 1 < mutStart_) return tpl_[i];

    // if we're beyond the mutation position, take the appropriate base
    else if (i > mutStart_)
        return tpl_[i - mutOff_];

    // otherwise if we're the base before mutStart_, 0, else 1 of our mutated tpl
    // params
    return mutTpl_[i == mutStart_];
}

bool Template::IsMutated() const { return mutated_; }

size_t VirtualTemplate::Length() const
{
    if (IsMutated()) return end_ - start_ + master_.mutOff_;

    return end_ - start_;
}

std::unique_ptr<AbstractRecursor> Template::CreateRecursor(std::unique_ptr<AbstractTemplate>&& tpl,
                                                           const MappedRead& mr,
                                                           double scoreDiff) const
{
    return cfg_->CreateRecursor(std::forward<std::unique_ptr<AbstractTemplate>>(tpl), mr,
                                scoreDiff);
}

double Template::SubstitutionRate(uint8_t prev, uint8_t curr) const
{
    return cfg_->SubstitutionRate(prev, curr);
}

const TemplatePosition& VirtualTemplate::operator[](const size_t i) const
{
    if (master_.IsMutated() && !pinStart_ && master_.mutEnd_ <= start_)
        return master_[start_ + i + master_.mutOff_];
    return master_[start_ + i];
}

bool VirtualTemplate::IsMutated() const
{
    return master_.IsMutated() && InRange(master_.mutStart_, master_.mutEnd_);
}

boost::optional<Mutation> VirtualTemplate::Mutate(const Mutation& mut)
{
    if (!master_.IsMutated()) throw std::runtime_error("virtual template badness");
    if (!InRange(mut.Start(), mut.End())) return boost::optional<Mutation>(boost::none);
    return boost::optional<Mutation>(Mutation(mut.Type, mut.Start() - start_, mut.Base));
}

std::unique_ptr<AbstractRecursor> VirtualTemplate::CreateRecursor(
    std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr, double scoreDiff) const
{
    return master_.CreateRecursor(std::forward<std::unique_ptr<AbstractTemplate>>(tpl), mr,
                                  scoreDiff);
}

double VirtualTemplate::SubstitutionRate(uint8_t prev, uint8_t curr) const
{
    return master_.SubstitutionRate(prev, curr);
}

}  // namespace Consensus
}  // namespace PacBio
