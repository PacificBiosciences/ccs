// Copyright (c) 2011-2014, Pacific Biosciences of California, Inc.
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

// Author: David Alexander

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include <pacbio/consensus/MonoMolecularIntegrator.h>
#include <pacbio/consensus/MultiMolecularIntegrator.h>
#include <pacbio/consensus/Mutation.h>

using std::string;
using std::vector;
using std::cout;
using std::endl;

namespace PacBio {
namespace Consensus {

vector<Mutation> Mutations(const AbstractIntegrator& ai, const size_t start, const size_t end);

vector<Mutation> Mutations(const AbstractIntegrator& ai);

vector<Mutation> NearbyMutations(vector<Mutation>* applied, vector<Mutation>* centers,
                                 const AbstractIntegrator& ai, const size_t neighborhood);
}
}

using namespace PacBio::Consensus;  // NOLINT

using ::testing::UnorderedElementsAreArray;

TEST(MutationEnumerationTest, TestAllMutationsMono)
{
    string tpl = "GAATC";
    MonoMolecularIntegrator ai(tpl, IntegratorConfig(), SNR(4, 4, 4, 4), "P6-C4");
    vector<Mutation> result = Mutations(ai);
    // 3 insertions, 3 substitutions, and 1 deletion per base
    //   and +4 for terminal insertions (1 beginning, 3 end)
    //   and -1 for the homopolymer AA deletion
    EXPECT_EQ(7 * tpl.length() + 4 - 1, result.size());
}

TEST(MutationEnumerationTest, TestAllMutationsMulti)
{
    string tpl = "GAATC";
    MultiMolecularIntegrator ai(tpl, IntegratorConfig());
    vector<Mutation> result = Mutations(ai);
    EXPECT_EQ(7 * tpl.length() + 4 - 1, result.size());
}

TEST(MutationEnumerationTest, TestNearbyMutations)
{
    string tpl = "GAATT";
    MonoMolecularIntegrator ai(tpl, IntegratorConfig(), SNR(4, 4, 4, 4), "P6-C4");

    vector<Mutation> centers = {Mutation(MutationType::SUBSTITUTION, 2, 'T')};
    vector<Mutation> result = NearbyMutations(&centers, &centers, ai, 1);
    // 7 for each of ATT,
    //   and +3 for terminal insertions (end)
    //   and -1 for hompolymer TT deletion
    EXPECT_EQ(7 * 3 + 3 - 1, result.size());

    centers.back() = Mutation(MutationType::SUBSTITUTION, 1, 'T');

    result = NearbyMutations(&centers, &centers, ai, 1);
    // 7 for each of ATT,
    //   and +4 for terminal insertions (1 beg, 3 end)
    //   and -1 for hompolymer TT deletion
    EXPECT_EQ(7 * 3 + 4 - 1, result.size());

    result = NearbyMutations(&centers, &centers, ai, 2);
    // 8 mutations for each of GAAT, except -2 for the homopolymer AA indels,
    //   and +3 for terminal insertions (not T, because it's not the first)
    EXPECT_EQ(7 * 4 + 4 - 1, result.size());

    centers.push_back(Mutation(MutationType::SUBSTITUTION, 3, 'G'));
    result = NearbyMutations(&centers, &centers, ai, 2);
    vector<Mutation> expected = Mutations(ai);
    EXPECT_EQ(expected.size(), result.size());
    EXPECT_THAT(result, UnorderedElementsAreArray(expected));
}
