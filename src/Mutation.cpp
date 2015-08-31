
#include <pacbio/consensus/Mutation.h>

namespace PacBio {
namespace Consensus {

Mutation::Mutation(MutationType type, size_t start, char base)
    : base_{base}
    , start_{start}
    , type_{type}
{ }

bool Mutation::IsDeletion() const
{
    return type_ == MutationType::DELETION;
}

bool Mutation::IsInsertion() const
{
    return type_ == MutationType::INSERTION;
}

bool Mutation::IsSubstitution() const
{
    return type_ == MutationType::SUBSTITUTION;
}

bool Mutation::IsAnyInsertion() const
{
    return type_ == MutationType::ANY_INSERTION;
}

bool Mutation::IsAnySubstitution() const
{
    return type_ == MutationType::ANY_SUBSTITUTION;
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
    if (type_ == MutationType::INSERTION ||
        type_ == MutationType::ANY_INSERTION)
        return start_;
    
    // if (type_ == MutationType::SUBSTITUTION ||
    //     type_ == MutationType::ANY_SUBSTITUTION ||
    //     type_ == MutationType::DELETION)
    return start_ + 1;
}

MutationType Mutation::Type() const
{
    return type_;
}

int Mutation::LengthDiff() const
{
    if (type_ == MutationType::SUBSTITUTION ||
        type_ == MutationType::ANY_SUBSTITUTION)
        return 0;

    if (type_ == MutationType::INSERTION ||
        type_ == MutationType::ANY_INSERTION)
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
