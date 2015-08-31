
#pragma once

#include <cassert>
#include <cmath>
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
    virtual const TemplatePosition& operator[](size_t i) const = 0;

    // virtual mutations (for mutation testing purposes)
    virtual bool IsMutated() const = 0;
    virtual void Mutate(const Mutation& m) = 0;
    virtual void Reset() = 0;

    // actually apply mutations
    virtual void ApplyMutation(const Mutation& mut) = 0;
    virtual void ApplyMutations(const std::vector<Mutation>& muts) = 0;

    // access model configuration
    virtual double BaseEmissionPr(char from, char to) const;
    virtual double CovEmissionPr(MoveType move, uint8_t cov) const;

    std::tuple<double, double> NormalParameters(size_t start, size_t end) const;

private:
    std::tuple<double, double> SiteNormalParameters(size_t i) const;
};

class Template : public AbstractTemplate
{
public:
    Template(const std::string& tpl, std::unique_ptr<ModelConfig>&& cfg)
        : cfg_(std::forward<std::unique_ptr<ModelConfig>>(cfg))
        , tpl_{cfg_->Populate(tpl)}
    { }

    const TemplatePosition& operator[](size_t i) const
    { return tpl_[i]; }

    bool IsMutated() const;
    void Mutate(const Mutation& mut);
    void Reset();

    void ApplyMutation(const Mutation& mut);
    void ApplyMutations(const std::vector<Mutation>& muts);

    inline
    double BaseEmissionPr(char from, char to) const
    {
        return cfg_->BaseEmissionPr(from, to);
    }

    inline
    double CovEmissionPr(MoveType move, uint8_t cov) const
    {
        return cfg_->CovEmissionPr(move, cov);        
    }

private:
    std::unique_ptr<ModelConfig> cfg_;
    std::vector<TemplatePosition> tpl_;

    friend class VirtualTemplate;
};

class VirtualTemplate : public AbstractTemplate
{
public:
    VirtualTemplate(const Template& master, size_t start, size_t end)
        : master_{&master}
        , start_{start}
        , end_{end}
    { assert(start_ < end_); }

    inline
    const TemplatePosition& operator[](size_t i) const
    { return (*master_)[i - start_]; }

    inline
    bool IsMutated() const
    {
        return master_->IsMutated();
    }

    inline
    void Mutate(const Mutation& m)
    {
        if (!master_->IsMutated())
            throw std::runtime_error("virtual template badness");
    }

    inline
    void Reset()
    {
        if (master_->IsMutated())
            throw std::runtime_error("virtual template badness");
    }

    inline
    void ApplyMutation(const Mutation& mut)
    { }

    inline
    void ApplyMutations(const std::vector<Mutation>& muts)
    { }

    inline
    double BaseEmissionPr(char from, char to) const
    {
        return master_->cfg_->BaseEmissionPr(from, to);
    }

    inline
    double CovEmissionPr(MoveType move, uint8_t cov) const
    {
        return master_->cfg_->CovEmissionPr(move, cov);
    }

private:
    Template const* master_;
    size_t start_;
    size_t end_;
};

} // namespace Consensus
} // namespace PacBio
