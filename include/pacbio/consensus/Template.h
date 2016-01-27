
#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <boost/optional.hpp>

#include <pacbio/consensus/ModelConfig.h>
#include <pacbio/consensus/Mutation.h>

namespace PacBio {
namespace Consensus {

class AbstractTemplate
{
public:
    virtual ~AbstractTemplate();

    virtual size_t Length() const = 0;
    virtual const TemplatePosition& operator[](size_t i) const = 0;

    // virtual mutations (for mutation testing purposes)
    virtual bool IsMutated() const = 0;
    virtual boost::optional<Mutation> Mutate(const Mutation& m) = 0;
    virtual void Reset() = 0;

    // actually apply mutations
    virtual void ApplyMutation(const Mutation& mut);
    virtual void ApplyMutations(std::vector<Mutation>* muts);

    // access model configuration
    virtual double BaseEmissionPr(MoveType move, char from, char to) const = 0;
    virtual double CovEmissionPr(MoveType move, uint8_t cov, char from, char to) const = 0;
    virtual double UndoCounterWeights(size_t nEmissions) const = 0;

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

    void ApplyMutation(const Mutation& mut);

    inline double BaseEmissionPr(MoveType move, char from, char to) const;
    inline double CovEmissionPr(MoveType move, uint8_t cov, char from, char to) const;
    inline double UndoCounterWeights(size_t nEmissions) const;

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
    void ApplyMutation(const Mutation& mut);

    inline double BaseEmissionPr(MoveType move, char from, char to) const;
    inline double CovEmissionPr(MoveType move, uint8_t cov, char from, char to) const;
    inline double UndoCounterWeights(size_t nEmissions) const;

private:
    Template const& master_;
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
double Template::BaseEmissionPr(MoveType move, char from, char to) const
{
    return cfg_->BaseEmissionPr(move, from, to);
}

double Template::CovEmissionPr(MoveType move, uint8_t cov, char from, char to) const
{
    return cfg_->CovEmissionPr(move, cov, from, to);
}

double Template::UndoCounterWeights(size_t nEmissions) const
{
    return cfg_->UndoCounterWeights(nEmissions);
}

size_t VirtualTemplate::Length() const
{
    if (IsMutated()) return end_ - start_ + master_.mutOff_;

    return end_ - start_;
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

double VirtualTemplate::BaseEmissionPr(MoveType move, char from, char to) const
{
    return master_.cfg_->BaseEmissionPr(move, from, to);
}

double VirtualTemplate::CovEmissionPr(MoveType move, uint8_t cov, char from, char to) const
{
    return master_.cfg_->CovEmissionPr(move, cov, from, to);
}

double VirtualTemplate::UndoCounterWeights(size_t nEmissions) const
{
    return master_.cfg_->UndoCounterWeights(nEmissions);
}

}  // namespace Consensus
}  // namespace PacBio
