// Author: David Alexander

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include <pacbio/consensus/Integrator.h>
#include <pacbio/consensus/Mutation.h>
#include <pacbio/consensus/Polish.h>
#include <pacbio/data/Read.h>

using std::string;
using std::vector;
using std::cout;
using std::endl;

using namespace PacBio::Data;

namespace PacBio {
namespace Consensus {

vector<Mutation> NearbyMutations(vector<Mutation>* applied, vector<Mutation>* centers,
                                 const Integrator& ai, const size_t neighborhood,
                                 const bool diploid = false);
}
}

using namespace PacBio::Consensus;  // NOLINT

using ::testing::UnorderedElementsAreArray;

TEST(MutationEnumerationTest, TestAllMutations)
{
    string tpl = "GAATC";
    Integrator ai(tpl, IntegratorConfig());
    vector<Mutation> result = Mutations(ai);
    EXPECT_EQ(7 * tpl.length() + 4 - 1, result.size());
}

TEST(MutationEnumerationTest, TestDiRepeatMutations)
{
    string tpl = "ACGTATATATACATATATTGCA";
    Integrator ai(tpl, IntegratorConfig());
    vector<Mutation> result = RepeatMutations(ai, RepeatConfig());
    ASSERT_EQ(4, result.size());
    EXPECT_EQ(Mutation::Insertion(3, "TA"), result[0]);
    EXPECT_EQ(Mutation::Deletion(3, 2), result[1]);
    EXPECT_EQ(Mutation::Insertion(12, "AT"), result[2]);
    EXPECT_EQ(Mutation::Deletion(12, 2), result[3]);
}

TEST(MutationEnumerationTest, TestTriRepeatMutations)
{
    string tpl = "ACGTCAGCAGCAGGAGGAGGTGCA";
    Integrator ai(tpl, IntegratorConfig());
    vector<Mutation> result = RepeatMutations(ai, RepeatConfig());
    ASSERT_EQ(4, result.size());
    EXPECT_EQ(Mutation::Insertion(4, "CAG"), result[0]);
    EXPECT_EQ(Mutation::Deletion(4, 3), result[1]);
    EXPECT_EQ(Mutation::Insertion(11, "AGG"), result[2]);
    EXPECT_EQ(Mutation::Deletion(11, 3), result[3]);
}

TEST(MutationEnumerationTest, TestNearbyMutations)
{
    string tpl = "GAATT";
    Integrator ai(tpl, IntegratorConfig());

    vector<Mutation> centers = {Mutation::Substitution(2, 'T')};
    vector<Mutation> result = NearbyMutations(&centers, &centers, ai, 1);
    // 7 for each of ATT,
    //   and +3 for terminal insertions (end)
    //   and -1 for hompolymer TT deletion
    EXPECT_EQ(7 * 3 + 3 - 1, result.size());

    centers.back() = Mutation::Substitution(1, 'T');

    result = NearbyMutations(&centers, &centers, ai, 1);
    // 7 for each of ATT,
    //   and +4 for terminal insertions (1 beg, 3 end)
    //   and -1 for hompolymer TT deletion
    EXPECT_EQ(7 * 3 + 4 - 1, result.size());

    result = NearbyMutations(&centers, &centers, ai, 2);
    // 8 mutations for each of GAAT, except -2 for the homopolymer AA indels,
    //   and +3 for terminal insertions (not T, because it's not the first)
    EXPECT_EQ(7 * 4 + 4 - 1, result.size());

    centers.push_back(Mutation::Substitution(3, 'G'));
    result = NearbyMutations(&centers, &centers, ai, 2);
    vector<Mutation> expected = Mutations(ai);
    EXPECT_EQ(expected.size(), result.size());
    EXPECT_THAT(result, UnorderedElementsAreArray(expected));
}
