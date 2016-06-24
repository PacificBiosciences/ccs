// Copyright (c) 2011-2016, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

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
