
#pragma once

#include <memory>
#include <utility>
#include <vector>

#include <pacbio/consensus/Read.h>
#include <pacbio/consensus/Template.h>

namespace PacBio {
namespace Consensus {

// forward declaration
class EvaluatorImpl;

class Evaluator
{
public:
    Evaluator(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr,
              double scoreDiff = 12.5);

    // move constructors
    Evaluator(Evaluator&&);

    ~Evaluator();

    size_t Length() const;
    char operator[](size_t i) const;
    operator std::string() const;

    double LL(const Mutation& mut);
    double LL() const;

    std::pair<double, double> NormalParameters() const;

    double ZScore() const;

    bool ApplyMutation(const Mutation& mut);
    bool ApplyMutations(std::vector<Mutation>* muts);

private:
    std::unique_ptr<EvaluatorImpl> impl_;
};

}  // namespace Consensus
}  // namespace PacBio
