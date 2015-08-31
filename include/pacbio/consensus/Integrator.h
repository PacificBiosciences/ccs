
#pragma once

#include <pacbio/consensus/Evaluator.h>
#include <pacbio/consensus/Exceptions.h>
#include <pacbio/consensus/Mutation.h>
#include <pacbio/consensus/ParameterTable.h>

namespace PacBio {
namespace Consensus {

class IntegratorConfig
{
public:
    IntegratorConfig(const ParameterTable& pt = ParameterTable::Default(),
                     double minZScore = -5.0);
    IntegratorConfig(double minZScore = -5.0);

    ParameterTable const* ParamTable;
    double MinZScore;
};

enum AddReadResult
{
    SUCCESS,
    ALPHA_BETA_MISMATCH,
    POOR_ZSCORE,
    OTHER
};

class AbstractIntegrator
{
public:
    virtual double LL(const Mutation& mut);
    
    double LL() const;
    double AvgZScore() const;

    virtual void ApplyMutation(const Mutation& mut) = 0;
    virtual void ApplyMutations(const std::vector<Mutation>& muts) = 0;

protected:
    AbstractIntegrator(const IntegratorConfig& cfg)
        : cfg_{cfg}
    { }

    AddReadResult AddRead(Evaluator&& eval); 
    
    IntegratorConfig cfg_;
    std::vector<Evaluator> evals_;
};

class MonoMolecularIntegrator : public AbstractIntegrator
{
public:
    MonoMolecularIntegrator(const std::string& tpl,
                            const IntegratorConfig& cfg,
                            const SNR& snr,
                            const std::string& model);

    AddReadResult AddRead(const MappedRead& read);

    double LL(const Mutation& mut);

private:
    std::string mdl_;
    Template tpl_;
};

class MultiMolecularIntegrator : public AbstractIntegrator
{
public:
    MultiMolecularIntegrator(const std::string& tpl,
                             const IntegratorConfig& cfg);

    AddReadResult AddRead(const MappedRead& read, const SNR& snr);

private:
    std::string tpl_;
};

} // namespace Consensus
} // namespace PacBio
