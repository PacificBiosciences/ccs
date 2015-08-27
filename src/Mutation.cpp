
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

Mutation::Mutation(MutationType type, size_t start, char base)
    : base_{base}
    , start_{start}
    , type_{type}
{ }

bool Mutation::IsDeletion() const
{
    return type_ == DELETION;
}

bool Mutation::IsInsertion() const
{
    return type_ == INSERTION;
}

bool Mutation::IsSubstitution() const
{
    return type_ == SUBSTITUTION;
}

bool Mutation::IsAnyInsertion() const
{
    return type_ == ANY_INSERTION;
}

bool Mutation::IsAnySubstitution() const
{
    return type_ == ANY_SUBSTITUTION;
}

char Mutation::Base() const
{
    return base_;
}

size_t Mutation::Start() const
{
    return start_;
}

size_t Mutation::End() const
{
    if (type_ == INSERTION || type_ == ANY_INSERTION)
        return start_;
    
    //if (type_ == SUBSTITUTION || type_ == ANY_SUBSTITUTION || type_ == DELETION)
    return start_ + 1;
}

MutationType Type() const
{
    return type_;
}

int Mutation::LengthDiff() const
{
    if (type_ == SUBSTITUTION || type_ == ANY_SUBSTITUTION)
        return 0;

    if (type_ == INSERTION || type_ == ANY_INSERTION)
        return 1;

    // type_ == deletion
    return -1;
}

ScoredMutation Mutation::WithScore(double score) const
{
    return ScoredMutation(*this, score);
}

ScoredMutation::ScoredMutation(const Mutation& mut, double score)
    : Mutation(mut)
    , score_{score}
{ }

double ScoredMutation::Score() const
{
    return score_;
}

} // namespace Consensus
} // namespace PacBio
