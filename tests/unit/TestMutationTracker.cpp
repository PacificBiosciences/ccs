// Author: David Seifert

#include <string>

#include <gtest/gtest.h>

#include <pacbio/consensus/Mutation.h>
#include "../src/MutationTracker.h"

using namespace PacBio::Consensus;

TEST(MutationTrackerTest, TestInterleavedMutations)
{
    // Test all combinations of different types of
    // mutations following each other. Furthermore
    // do it in two rounds, such that everything is
    // interleaved in the most complex way.

    // short template from all4mers
    //
    // 0    5    10   15   20
    // ATAATCAGCGACCTCCTAGCCAGTC
    MutationTracker mutTestTracker{"ATAATCAGCGACCTCCTAGCCAGTC"};

    // 1st round of Mutations
    //
    //  Original:  ATA ATCAGC GACCTCCTAGCCAGTC
    //              S I      I  X     X S
    // 1st Round:  AGAcATCAGCtGA-CTCCT-GTCAGTC
    //              | |      |  |     | |
    //            S,1 |    I,9  |  D,17 |
    //               I,3      D,11     S,19
    //
    std::vector<Mutation> firstRoundMutations{
        Mutation::Substitution(1, 'G'), Mutation::Insertion(3, 'C'),
        Mutation::Insertion(9, 'T'),    Mutation::Deletion(11, 1),
        Mutation::Deletion(17, 1),      Mutation::Substitution(19, 'T')};
    std::sort(firstRoundMutations.begin(), firstRoundMutations.end(), Mutation::SiteComparer);
    mutTestTracker.AddSortedMutations(firstRoundMutations);

    // 2nd round of Mutations
    //
    //  Original:  ATA AT CAGC GACCTC CTAGCCAGTC
    //              s i       i  x      x s
    // 1st Round:  AGAcAT CAGCtGA-CTC CT-GTCAGTC
    //                   I  S      X I      S X
    // 2nd Round:  AGAcATgCATCtGA-C-CgCT-GTCTG-C
    //                   |  |      | |      | |
    //                 I,6  |   D,14 |   S,21 |
    //                     S,8      I,16     D,23
    //
    std::vector<Mutation> secondRoundMutations{
        Mutation::Insertion(6, 'G'),  Mutation::Substitution(8, 'T'),  Mutation::Deletion(14, 1),
        Mutation::Insertion(16, 'G'), Mutation::Substitution(21, 'T'), Mutation::Deletion(23, 1)};
    std::sort(secondRoundMutations.begin(), secondRoundMutations.end(), Mutation::SiteComparer);
    mutTestTracker.AddSortedMutations(secondRoundMutations);

    const auto finalMapping = mutTestTracker.MappingToOriginalTpl();

#if 0
    // If you ever need to debug...
    constexpr const std::array<const char*, 3> FullMutationTypes = {{"DEL", "INS", "SUB"}};
    for (const auto& i : finalMapping) {
        std::cout << FullMutationTypes[static_cast<uint8_t>(i.mutType)] << '\t' << i.pos;
        for (const auto& i : i.mutants) {
            std::cout << '\t' << i;
        }
        std::cout << std::endl;
    }
#endif

    const std::vector<DiploidSite> correctOutput = {
        {MutationType::SUBSTITUTION, std::vector<char>{'G'}, 1},
        {MutationType::INSERTION, std::vector<char>{'C'}, 3},
        {MutationType::INSERTION, std::vector<char>{'G'}, 5},
        {MutationType::SUBSTITUTION, std::vector<char>{'T'}, 7},
        {MutationType::INSERTION, std::vector<char>{'T'}, 9},
        {MutationType::DELETION, std::vector<char>{}, 11},
        {MutationType::DELETION, std::vector<char>{}, 13},
        {MutationType::INSERTION, std::vector<char>{'G'}, 15},
        {MutationType::DELETION, std::vector<char>{}, 17},
        {MutationType::SUBSTITUTION, std::vector<char>{'T'}, 19},
        {MutationType::SUBSTITUTION, std::vector<char>{'T'}, 21},
        {MutationType::DELETION, std::vector<char>{}, 23}};

    EXPECT_EQ(correctOutput, finalMapping);
}

TEST(MutationTrackerTest, TestFrontDeletion)
{
    // Test that we catch deletions at the beginning
    MutationTracker mutTestTracker{"AACCGGTT"};

    std::vector<Mutation> Mutations{Mutation::Deletion(0, 2)};
    std::sort(Mutations.begin(), Mutations.end(), Mutation::SiteComparer);
    mutTestTracker.AddSortedMutations(Mutations);

    const auto finalMapping = mutTestTracker.MappingToOriginalTpl();

    const std::vector<DiploidSite> correctOutput = {
        {MutationType::DELETION, std::vector<char>{}, 0},
        {MutationType::DELETION, std::vector<char>{}, 1}};

    EXPECT_EQ(correctOutput, finalMapping);
}

TEST(MutationTrackerTest, TestBackDeletion)
{
    // Test that we catch deletions at the end
    MutationTracker mutTestTracker{"AACCGGTT"};

    std::vector<Mutation> Mutations{Mutation::Deletion(6, 2)};
    std::sort(Mutations.begin(), Mutations.end(), Mutation::SiteComparer);
    mutTestTracker.AddSortedMutations(Mutations);

    const auto finalMapping = mutTestTracker.MappingToOriginalTpl();

    const std::vector<DiploidSite> correctOutput = {
        {MutationType::DELETION, std::vector<char>{}, 6},
        {MutationType::DELETION, std::vector<char>{}, 7}};

    EXPECT_EQ(correctOutput, finalMapping);
}

TEST(MutationTrackerTest, TestInsertionSubstitution)
{
    // Test that we catch deletions at the end
    MutationTracker mutTestTracker{"AT"};

    std::vector<Mutation> firstRoundMutations{Mutation::Insertion(1, "GG")};
    std::sort(firstRoundMutations.begin(), firstRoundMutations.end(), Mutation::SiteComparer);
    mutTestTracker.AddSortedMutations(firstRoundMutations);

    std::vector<Mutation> secondRoundMutations{Mutation::Substitution(1, "CC")};
    std::sort(secondRoundMutations.begin(), secondRoundMutations.end(), Mutation::SiteComparer);
    mutTestTracker.AddSortedMutations(secondRoundMutations);

    const auto finalMapping = mutTestTracker.MappingToOriginalTpl();

    const std::vector<DiploidSite> correctOutput = {
        {MutationType::INSERTION, std::vector<char>{'C'}, 1},
        {MutationType::INSERTION, std::vector<char>{'C'}, 1}};

    EXPECT_EQ(correctOutput, finalMapping);
}
