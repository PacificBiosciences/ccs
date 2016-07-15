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
