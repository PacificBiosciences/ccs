
#pragma once

#include <cstdint>
#include <functional>
#include <iostream>
#include <memory>
#include <set>

#include <pacbio/consensus/Evaluator.h>
#include <pacbio/consensus/Exceptions.h>
#include <pacbio/consensus/Mutation.h>

namespace PacBio {
namespace Consensus {

struct IntegratorConfig
{
    double MinZScore;
    double ScoreDiff;

    IntegratorConfig(double minZScore = -3.5, double scoreDiff = 12.5);
};

enum struct AddReadResult : uint8_t
{
    SUCCESS,
    ALPHA_BETA_MISMATCH,
    POOR_ZSCORE,
    OTHER,
    /*
     This enum is used in other places to
     both size an array and index into it.  So
     users can know the number of elements needed in
     an array that could count valus of this enum, this
     should always be the last value in the enum.
     */
    SIZE
};

inline std::ostream& operator<<(std::ostream& os, AddReadResult result)
{
    static const char* names[] = {"SUCCESS", "ALPHA/BETA MISMATCH", "POOR Z-SCORE", "OTHER"};
    os << names[static_cast<size_t>(result)];
    return os;
}

std::set<std::string> SupportedChemistries();

class AbstractIntegrator
{
public:
    virtual ~AbstractIntegrator();

    virtual size_t TemplateLength() const = 0;

    virtual char operator[](size_t i) const = 0;
    virtual operator std::string() const = 0;

    virtual double LL(const Mutation& mut);
    virtual double LL() const;

    double AvgZScore() const;
    std::vector<double> ZScores() const;
    std::vector<std::pair<double, double>> NormalParameters() const;

    virtual void ApplyMutation(const Mutation& mut) = 0;
    virtual void ApplyMutations(std::vector<Mutation>* muts) = 0;

    virtual AddReadResult AddRead(const MappedRead& read) = 0;

    // For debugging purposes
    // (Note that these include results include all evaluators, even the inactive ones)
    std::vector<double> LLs(const Mutation& mut);
    std::vector<double> LLs() const;
    std::vector<std::string> ReadNames() const;

protected:
    Mutation ReverseComplement(const Mutation& mut) const;

    AbstractIntegrator(const IntegratorConfig& cfg);

    // move constructor
    AbstractIntegrator(AbstractIntegrator&&);

    AddReadResult AddRead(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& read);

    IntegratorConfig cfg_;
    std::vector<Evaluator> evals_;
};

class MonoMolecularIntegrator : public AbstractIntegrator
{
public:
    MonoMolecularIntegrator(const std::string& tpl, const IntegratorConfig& cfg, const SNR& snr,
                            const std::string& model);

    // move constructor
    MonoMolecularIntegrator(MonoMolecularIntegrator&&);

    size_t TemplateLength() const;

    char operator[](size_t i) const;
    operator std::string() const;

    double LL(const Mutation& mut);
    inline double LL() const { return AbstractIntegrator::LL(); }
    void ApplyMutation(const Mutation& mut);
    void ApplyMutations(std::vector<Mutation>* muts);

    AddReadResult AddRead(const MappedRead& read);

protected:
    std::string mdl_;
    SNR snr_;
    Template fwdTpl_;
    Template revTpl_;
};

class MultiMolecularIntegrator : public AbstractIntegrator
{
public:
    MultiMolecularIntegrator(const std::string& tpl, const IntegratorConfig& cfg);

    size_t TemplateLength() const;

    char operator[](size_t i) const;
    operator std::string() const;

    void ApplyMutation(const Mutation& mut);
    void ApplyMutations(std::vector<Mutation>* muts);

    AddReadResult AddRead(const MappedRead& read);

protected:
    std::unique_ptr<AbstractTemplate> GetTemplate(const MappedRead& read, const SNR& snr);

    std::string fwdTpl_;
    std::string revTpl_;

private:
    friend struct std::hash<MultiMolecularIntegrator>;
};

}  // namespace Consensus
}  // namespace PacBio
