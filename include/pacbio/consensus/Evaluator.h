
#pragma once

#include <memory>
#include <utility>
#include <vector>

#include <pacbio/consensus/Read.h>
#include <pacbio/consensus/Template.h>

namespace PacBio {
namespace Consensus {

enum struct EvaluatorState : uint8_t
{
    VALID,
    ALPHA_BETA_MISMATCH,
    POOR_ZSCORE,
    NULL_TEMPLATE,
    DISABLED
};

// forward declaration
class EvaluatorImpl;

class Evaluator
{
public:
    Evaluator(EvaluatorState);
    Evaluator(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr, double minZScore,
              double scoreDiff);

    // move constructors
    Evaluator(Evaluator&&);
    Evaluator& operator=(Evaluator&&);

    ~Evaluator();

    size_t Length() const;
    StrandEnum Strand() const;

    operator bool() const;
    operator std::string() const;

    double LL(const Mutation& mut);
    double LL() const;

    std::pair<double, double> NormalParameters() const;

    double ZScore() const;

    bool ApplyMutation(const Mutation& mut);
    bool ApplyMutations(std::vector<Mutation>* muts);

    EvaluatorState Status() const;
    void Release();

private:
    void CheckInvariants();

private:
    std::unique_ptr<EvaluatorImpl> impl_;
    EvaluatorState state_;
};

}  // namespace Consensus
}  // namespace PacBio
