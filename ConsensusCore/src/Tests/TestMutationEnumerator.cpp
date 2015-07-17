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
#include <boost/assign.hpp>
#include <boost/assign/std/set.hpp>
#include <boost/range/as_array.hpp>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include <ConsensusCore/Utils.hpp>
#include <ConsensusCore/Mutation.hpp>
#include <ConsensusCore/MutationEnumerator.hpp>

using std::string;
using std::vector;
using std::cout;
using std::endl;

using namespace boost::assign;  // NOLINT
using namespace ConsensusCore;  // NOLINT

using ::testing::UnorderedElementsAreArray;

TEST(MutationEnumerationTest, TestAllMutations)
{
    std::string tpl = "GAATC";
    std::vector<Mutation> result = AllSingleBaseMutationEnumerator(tpl).Mutations();
    // 4 insertions, 3 substitutions, and 1 deletion per base
    EXPECT_EQ(8*tpl.length(), result.size());
}

TEST(MutationEnumerationTest, TestUniqueMutations)
{
    std::string tpl = "GAATC";
    std::vector<Mutation> result = UniqueSingleBaseMutationEnumerator(tpl).Mutations();
    // 3 insertions, 3 substitions, and 1 deletion per base,
    // except the first (which has an extra insertion),
    // and the homopolymeric A (which is less a deletion)
    EXPECT_EQ(7*tpl.length() + 1 - 1, result.size());
}


TEST(MutationEnumerationTest, TestUniqueNearbyMutations)
{
    std::string tpl = "GAATC";

    std::vector<Mutation> centers;
    centers.push_back(Mutation(SUBSTITUTION, 1, 'T'));

    UniqueSingleBaseMutationEnumerator enumerator(tpl);
    std::vector<Mutation> result = UniqueNearbyMutations(enumerator, centers, 1);
    // 8 mutations for the G,
    // but only 7 for the A because we don't want a repeat insertion
    EXPECT_EQ(8 + 7, result.size());

    result = UniqueNearbyMutations(enumerator, centers, 2);
    // 8 for the first, 7 for the second, 6 for the third (no homopolymeric deletion)
    EXPECT_EQ(8 + 7 + 6, result.size());

    centers.push_back(Mutation(SUBSTITUTION, 3, 'G'));
    result = UniqueNearbyMutations(enumerator, centers, 2);
    std::vector<Mutation> expected = UniqueSingleBaseMutationEnumerator(tpl).Mutations();
    EXPECT_THAT(result, UnorderedElementsAreArray(expected));
}


TEST(MutationEnumerationTest, TestDinucleotideMutations)
{
    std::string tpl = "ACACACGCGCGTGTG";
    std::vector<Mutation> result = DinucleotideRepeatMutationEnumerator(tpl, 3).Mutations();
    // 4 extra mutations because of ACACAC, CGCGCG, but not GTGTG
    EXPECT_EQ(4, result.size());

    std::vector<Mutation> expected;
    expected.push_back(Mutation(INSERTION, 0, 0, std::string("AC")));
    expected.push_back(Mutation(DELETION, 0, 2, std::string("")));
    expected.push_back(Mutation(INSERTION, 5, 5, std::string("CG")));
    expected.push_back(Mutation(DELETION, 5, 7, std::string("")));
    EXPECT_THAT(result, UnorderedElementsAreArray(expected));
}


TEST(MutationEnumerationTest, TestTrinucleotideMutations)
{

    std::string tpl = "ACAACAACAGCAGCAGTAGTAG";
    std::vector<Mutation> result = RepeatMutationEnumerator(tpl, 3, 3).Mutations();
    // 4 extra mutations because of ACAACAACA, CAGCAGCAG, but not AGTAGTAG
    EXPECT_EQ(4, result.size());

    std::vector<Mutation> expected;
    expected.push_back(Mutation(INSERTION, 0, 0, std::string("ACA")));
    expected.push_back(Mutation(DELETION, 0, 3, std::string("")));
    expected.push_back(Mutation(INSERTION, 7, 7, std::string("CAG")));
    expected.push_back(Mutation(DELETION, 7, 10, std::string("")));
    EXPECT_THAT(result, UnorderedElementsAreArray(expected));
}
