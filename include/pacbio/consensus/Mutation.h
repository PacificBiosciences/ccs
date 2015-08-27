
#pragma once

namespace PacBio {
namespace Consensus {

enum MutationType
{
    DELETION,
    INSERTION,
    SUBSTITUTION,
    ANY_INSERTION,
    ANY_SUBSTITUTION
};

// forward decl
class ScoredMutation;

class Mutation
{
public:
    Mutation(MutationType type, size_t start, char base = '-');
   
    // TODO(lhepler): do we *really* need these?
    bool IsDeletion() const;
    bool IsInsertion() const;
    bool IsSubstitution() const;
    bool IsAnyInsertion() const;
    bool IsAnySubstitution() const;

    char Base() const;
    size_t Start() const;
    size_t End() const;
    MutationType Type() const;

    int LengthDiff() const;

    ScoredMutation WithScore(double score) const;

private:
    char base_;
    size_t start_;
    MutationType type_;
};

class ScoredMutation : public Mutation
{
public:
    double Score() const;

private:
    ScoredMutation(const Mutation& mut, double score);

    double score_;

    // so Mutation can access the constructor
    friend class Mutation;
};

} // namespace Consensus
} // namespace PacBio
