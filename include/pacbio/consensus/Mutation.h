// Author: Lance Hepler

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <ostream>
#include <string>
#include <tuple>
#include <vector>

#include <boost/numeric/conversion/cast.hpp>
#include <boost/optional.hpp>

// Initialize data structures, do NOT remove
#include <pacbio/consensus/internal/ModelInternalInitializer.h>

namespace PacBio {
namespace Consensus {

// fwd decls
class ScoredMutation;

enum struct MutationType : uint8_t
{
    DELETION,
    INSERTION,
    SUBSTITUTION
};

/// the region position and length, and any alternative bases,
/// and the score too if it comes to that
class Mutation
{
public:
    Mutation(const Mutation&) = default;
    Mutation(Mutation&&) = default;

    Mutation& operator=(const Mutation&) = default;
    Mutation& operator=(Mutation&&) = default;

    // named constructors
    static Mutation Deletion(size_t start, size_t length);
    static Mutation Insertion(size_t start, char base);
    static Mutation Insertion(size_t start, std::string bases);
    static Mutation Substitution(size_t start, char base);
    static Mutation Substitution(size_t start, std::string bases);

    bool IsDeletion() const { return Type() == MutationType::DELETION; };
    bool IsInsertion() const { return Type() == MutationType::INSERTION; }
    bool IsSubstitution() const { return Type() == MutationType::SUBSTITUTION; }

    virtual bool IsScored() const { return false; }

    size_t Start() const { return start_; }
    size_t End() const { return start_ + length_; }
    size_t Length() const { return length_; }

    /// Returns the length difference introduced by this mutation.
    int LengthDiff() const
    {
        return boost::numeric_cast<int>(Bases().size()) - boost::numeric_cast<int>(length_);
    };

    size_t EditDistance() const { return std::max(Bases().size(), length_); }

    const std::string& Bases() const { return bases_; }
    MutationType Type() const { return type_; }

    /// Projects the Mutation from its original coordinates into the coordinate region
    /// defined by (start, length). If the mutation is outside this region, return none.
    boost::optional<Mutation> Translate(size_t start, size_t length) const;

    virtual bool operator==(const Mutation& other) const
    {
        if (Type() != other.Type() || Start() != other.Start()) return false;
        if (IsDeletion()) return Length() == other.Length();
        // insertion or substitution
        return Bases() == other.Bases();
    }

    virtual operator std::string() const;

    /// Uses this and the provided score to create and return a ScoredMutation.
    ScoredMutation WithScore(double score) const;

    /// Comparer to sort mutations by start/end.
    static bool SiteComparer(const Mutation& lhs, const Mutation& rhs)
    {
        // perform a lexicographic sort on End, Start, IsDeletion
        //   Deletions override everybody, so they get applied last,
        //   which implies, because false < true, use of !IsDeletion
        const auto l = std::make_tuple(lhs.End(), lhs.Start(), !lhs.IsDeletion());
        const auto r = std::make_tuple(rhs.End(), rhs.Start(), !rhs.IsDeletion());
        return l < r;
    }

    /// Used by the diploid API to store pvalues
    inline Mutation WithPvalue(const double pvalue)
    {
        auto result = *this;
        result.pvalue_ = pvalue;
        return result;
    }

    inline boost::optional<double> GetPvalue() const { return pvalue_; }

private:
    std::string bases_;
    MutationType type_;
    size_t start_;
    size_t length_;
    boost::optional<double> pvalue_;

    Mutation(MutationType type, size_t start, size_t length);
    Mutation(MutationType type, size_t start, std::string bases);
    Mutation(MutationType type, size_t start, char);
};

class ScoredMutation : public Mutation
{
public:
    double Score;

    ScoredMutation(const ScoredMutation&) = default;
    ScoredMutation(ScoredMutation&&) = default;

    ScoredMutation& operator=(const ScoredMutation&) = default;
    ScoredMutation& operator=(ScoredMutation&&) = default;

    virtual bool IsScored() const override { return true; }
    virtual bool operator==(const ScoredMutation& other) const
    {
        return Score == other.Score && static_cast<const Mutation&>(*this) == other;
    }
    virtual bool operator==(const Mutation& other) const override { return false; }

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
}
}  // ::PacBio::Consensus
