
#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

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
    virtual void Mutate(const Mutation& m) = 0;
    virtual void Reset() = 0;

    // actually apply mutations
    virtual void ApplyMutation(const Mutation& mut) = 0;
    virtual void ApplyMutations(std::vector<Mutation>* muts);

    // access model configuration
    virtual double BaseEmissionPr(MoveType move, char from, char to) const = 0;
    virtual double CovEmissionPr(MoveType move, uint8_t cov) const = 0;
    virtual double UndoCounterWeights(size_t nEmissions) const = 0;

    std::tuple<double, double> NormalParameters(size_t start, size_t end) const;

private:
    std::tuple<double, double> SiteNormalParameters(size_t i) const;
};

class Template : public AbstractTemplate
{
public:
    Template(const std::string& tpl, std::unique_ptr<ModelConfig>&& cfg);

    size_t Length() const;
    const TemplatePosition& operator[](size_t i) const;

    bool IsMutated() const;
    void Mutate(const Mutation& mut);
    void Reset();

    void ApplyMutation(const Mutation& mut);

    inline double BaseEmissionPr(MoveType move, char from, char to) const
    {
        return cfg_->BaseEmissionPr(move, from, to);
    }

    inline double CovEmissionPr(MoveType move, uint8_t cov) const
    {
        return cfg_->CovEmissionPr(move, cov);
    }

    inline double UndoCounterWeights(size_t nEmissions) const
    {
        return cfg_->UndoCounterWeights(nEmissions);
    }

private:
    std::unique_ptr<ModelConfig> cfg_;
    std::vector<TemplatePosition> tpl_;
    bool mutated_;
    size_t mutPos_;
    int mutOff_;
    TemplatePosition mutTpl_[2];

    friend class VirtualTemplate;
};

class VirtualTemplate : public AbstractTemplate
{
public:
    VirtualTemplate(const Template& master, size_t start, size_t end);

    inline size_t Length() const { return end_ - start_ + master_.mutOff_; }
    inline const TemplatePosition& operator[](size_t i) const
    {
        return master_[i - start_];
    }

    inline bool IsMutated() const { return master_.IsMutated(); }
    inline void Mutate(const Mutation& m)
    {
        if (!master_.IsMutated())
            throw std::runtime_error("virtual template badness");
    }

    inline void Reset() {}
    // I would like to do the following, but I don't reset the master
    //   before resetting the children
    // { if (master_->IsMutated()) throw std::runtime_error("virtual template
    // badness"); }

    void ApplyMutation(const Mutation& mut);

    inline double BaseEmissionPr(MoveType move, char from, char to) const
    {
        return master_.cfg_->BaseEmissionPr(move, from, to);
    }

    inline double CovEmissionPr(MoveType move, uint8_t cov) const
    {
        return master_.cfg_->CovEmissionPr(move, cov);
    }

    inline double UndoCounterWeights(size_t nEmissions) const
    {
        return master_.cfg_->UndoCounterWeights(nEmissions);
    }

private:
    Template const& master_;
    size_t start_;
    size_t end_;
};

}  // namespace Consensus
}  // namespace PacBio
