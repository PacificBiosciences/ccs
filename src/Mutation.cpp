// Author: Lance Hepler

#include <algorithm>
#include <cassert>
#include <sstream>
#include <stdexcept>
#include <utility>

#include <pacbio/consensus/Mutation.h>

namespace PacBio {
namespace Consensus {

Mutation Mutation::Deletion(size_t start, size_t length)
{
    return Mutation(MutationType::DELETION, start, length);
}

Mutation Mutation::Insertion(size_t start, const char base)
{
    return Mutation(MutationType::INSERTION, start, base);
}

Mutation Mutation::Insertion(size_t start, std::string bases)
{
    return Mutation(MutationType::INSERTION, start, std::move(bases));
}

Mutation Mutation::Substitution(size_t start, const char base)
{
    return Mutation(MutationType::SUBSTITUTION, start, base);
}

Mutation Mutation::Substitution(size_t start, std::string bases)
{
    return Mutation(MutationType::SUBSTITUTION, start, std::move(bases));
}

boost::optional<Mutation> Mutation::Translate(size_t start, size_t length) const
{
    // if the mutation end is before our start, or our end
    //   is less than the mutation start, we don't overlap the mutation:
    //   template:       [---)
    //   mutation:   [---)
    //   mutation:           [---)
    if (End() + IsInsertion() < start || (start + length + IsInsertion()) <= Start())
        return boost::none;
    // what remains is now one of three/five possibilities each:
    //   template:     [-------)
    //   mutation:   [---)-------)
    //   mutation:       [---)---)
    //   mutation:           [---)
    //   start = max ^ ^
    //   end =           min ^ ^
    const size_t newStart = std::max(Start(), start);
    const size_t newLen = std::min(End(), start + length) - newStart;
    if (IsInsertion()) return Mutation::Insertion(newStart - start, Bases());
    if (newLen == 0) return boost::none;
    if (IsDeletion()) return Mutation::Deletion(newStart - start, newLen);
    return Mutation::Substitution(newStart - start, Bases().substr(newStart - Start(), newLen));
}

Mutation::operator std::string() const
{
    std::ostringstream ss;
    ss << (*this);
    return ss.str();
}

ScoredMutation Mutation::WithScore(double score) const { return ScoredMutation(*this, score); }

Mutation::Mutation(const MutationType type, const size_t start, const size_t length)
    : bases_(), type_{type}, start_{start}, length_{length}
{
    assert(type == MutationType::DELETION);
}

Mutation::Mutation(const MutationType type, const size_t start, const char base)
    : bases_(1, base)
    , type_{type}
    , start_{start}
    , length_{(type_ == MutationType::INSERTION) ? size_t(0) : size_t(1)}
{
}

Mutation::Mutation(const MutationType type, const size_t start, std::string bases)
    : bases_{std::move(bases)}
    , type_{type}
    , start_{start}
    , length_{(type_ == MutationType::INSERTION) ? size_t(0) : bases_.length()}
{
    assert(bases_.length() > 0);
}

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
        default:
            throw std::invalid_argument("invalid MutationType");
    }
    return out;
}

std::ostream& operator<<(std::ostream& out, const Mutation& mut)
{
    if (mut.IsDeletion())
        return out << "Mutation::Deletion(" << mut.Start() << ", " << mut.Length() << ')';
    else if (mut.IsInsertion())
        return out << "Mutation::Insertion(" << mut.Start() << ", \"" << mut.Bases() << "\")";
    // mut.IsSubstitution()
    return out << "Mutation::Substitution(" << mut.Start() << ", \"" << mut.Bases() << "\")";
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
        if (it->IsInsertion())
            newTpl.insert(it->Start(), it->Bases());
        else
            newTpl.replace(it->Start(), it->Length(), it->Bases());
    }

    return newTpl;
}

}  // namespace Consensus
}  // namespace PacBio
