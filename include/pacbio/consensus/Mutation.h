
#pragma once

#include <cstddef>
#include <cstdint>
#include <ostream>
#include <string>
#include <tuple>
#include <vector>

namespace PacBio {
namespace Consensus {

enum struct MutationType : uint8_t
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
    char Base;
    MutationType Type;

    Mutation(MutationType type, size_t start, char base = '-');

    // TODO(lhepler): do we *really* need these?
    bool IsDeletion() const;
    bool IsInsertion() const;
    bool IsSubstitution() const;
    bool IsAnyInsertion() const;
    bool IsAnySubstitution() const;

    size_t Start() const;
    size_t End() const;

    int LengthDiff() const;

    bool operator==(const Mutation& other) const;
    operator std::string() const;

    ScoredMutation WithScore(double score) const;

    static bool SiteComparer(const Mutation& lhs, const Mutation& rhs)
    {
        // perform a lexicographic sort on End, Start, IsDeletion
        const auto l = std::make_tuple(lhs.End(), lhs.Start(), lhs.IsDeletion());
        const auto r = std::make_tuple(rhs.End(), rhs.Start(), rhs.IsDeletion());
        return l < r;
    }

private:
    size_t start_;
};

class ScoredMutation : public Mutation
{
public:
    double Score;

    static bool ScoreComparer(const ScoredMutation& lhs, const ScoredMutation& rhs)
    {
        return lhs.Score < rhs.Score;
    }

private:
    ScoredMutation(const Mutation& mut, double score);

    // so Mutation can access the constructor
    friend class Mutation;
};

std::ostream& operator<<(std::ostream& out, MutationType type);
std::ostream& operator<<(std::ostream& out, const Mutation& mut);
std::ostream& operator<<(std::ostream& out, const ScoredMutation& smut);

std::string ApplyMutations(const std::string& tpl, std::vector<Mutation>* muts);

}  // namespace Consensus
}  // namespace PacBio
