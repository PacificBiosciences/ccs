
#include <algorithm>
#include <cassert>
#include <sstream>
#include <stdexcept>

#include <pacbio/consensus/Mutation.h>

namespace PacBio {
namespace Consensus {

Mutation::Mutation(MutationType type, size_t start, char base)
    : Base{base}, Type{type}, start_{start}
{
    assert(Type == MutationType::DELETION || Type == MutationType::ANY_INSERTION ||
           Type == MutationType::ANY_SUBSTITUTION || Base == 'A' || Base == 'C' || Base == 'G' ||
           Base == 'T');
}

bool Mutation::IsDeletion() const { return Type == MutationType::DELETION; }
bool Mutation::IsInsertion() const { return Type == MutationType::INSERTION; }
bool Mutation::IsSubstitution() const { return Type == MutationType::SUBSTITUTION; }
bool Mutation::IsAnyInsertion() const { return Type == MutationType::ANY_INSERTION; }
bool Mutation::IsAnySubstitution() const { return Type == MutationType::ANY_SUBSTITUTION; }
size_t Mutation::Start() const { return start_; }
size_t Mutation::End() const
{
    if (Type == MutationType::INSERTION || Type == MutationType::ANY_INSERTION) return start_;

    // if (Type == MutationType::SUBSTITUTION ||
    //     Type == MutationType::ANY_SUBSTITUTION ||
    //     Type == MutationType::DELETION)
    return start_ + 1;
}

int Mutation::LengthDiff() const
{
    if (Type == MutationType::SUBSTITUTION || Type == MutationType::ANY_SUBSTITUTION) return 0;

    if (Type == MutationType::INSERTION || Type == MutationType::ANY_INSERTION) return 1;

    // Type == deletion
    return -1;
}

bool Mutation::operator==(const Mutation& other) const
{
    return Type == other.Type && Base == other.Base && start_ == other.start_;
}

Mutation::operator std::string() const
{
    std::stringstream ss;
    ss << (*this);
    return ss.str();
}

ScoredMutation Mutation::WithScore(double score) const { return ScoredMutation(*this, score); }
ScoredMutation::ScoredMutation(const Mutation& mut, double score) : Mutation(mut), Score{score} {}
std::ostream& operator<<(std::ostream& out, const MutationType type)
{
    out << "MutationType::";
    switch (type) {
        case MutationType::DELETION:
            out << "DELETION";
            break;
        case MutationType::INSERTION:
            out << "INSERTION";
            break;
        case MutationType::SUBSTITUTION:
            out << "SUBSTITUTION";
            break;
        case MutationType::ANY_INSERTION:
            out << "ANY_INSERTION";
            break;
        case MutationType::ANY_SUBSTITUTION:
            out << "ANY_SUBSTITUTION";
            break;
        default:
            throw std::invalid_argument("invalid MutationType");
    }
    return out;
}

std::ostream& operator<<(std::ostream& out, const Mutation& mut)
{
    return out << "Mutation(" << mut.Type << ", " << mut.Start() << ", '" << mut.Base << "')";
}

std::ostream& operator<<(std::ostream& out, const ScoredMutation& smut)
{
    return out << "ScoredMutation(" << static_cast<Mutation>(smut) << ", '" << smut.Score << "')";
}

std::string ApplyMutations(const std::string& oldTpl, std::vector<Mutation>* const muts)
{
    std::sort(muts->begin(), muts->end(), Mutation::SiteComparer);
    std::vector<Mutation>::const_reverse_iterator it;

    if (muts->empty() || oldTpl.empty()) return oldTpl;

    // TODO(lhepler) make this algorithm not (n^2)
    std::string newTpl(oldTpl);

    for (it = muts->crbegin(); it != muts->crend(); ++it) {
        if (it->IsDeletion())
            newTpl.erase(newTpl.begin() + it->Start());
        else if (it->IsInsertion())
            newTpl.insert(newTpl.begin() + it->Start(), it->Base);
        else if (it->IsSubstitution())
            newTpl[it->Start()] = it->Base;
    }

    return newTpl;
}

}  // namespace Consensus
}  // namespace PacBio
