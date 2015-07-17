// Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.
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

#include <gtest/gtest.h>
#include <boost/assign.hpp>
#include <string>
#include <vector>

#include <ConsensusCore/Features.hpp>
#include <ConsensusCore/Mutation.hpp>
#include <ConsensusCore/Read.hpp>
#include <ConsensusCore/Sequence.hpp>
#include <ConsensusCore/Quiver/MutationScorer.hpp>
#include <ConsensusCore/Quiver/QuiverConfig.hpp>
#include <ConsensusCore/Quiver/QvEvaluator.hpp>
#include <ConsensusCore/Quiver/ReadScorer.hpp>
#include <ConsensusCore/Quiver/SimpleRecursor.hpp>
#include <ConsensusCore/Quiver/SseRecursor.hpp>

#include "ParameterSettings.hpp"

using namespace ConsensusCore;  // NOLINT
using namespace boost::assign;  // NOLINT

typedef testing::Types<SimpleQvRecursor,
                       SseQvRecursor,
                       SparseSimpleQvRecursor,
                       SparseSseQvRecursor>    AllRecursorTypes;

TYPED_TEST_CASE(MutationScorerTest,          AllRecursorTypes);


//
// ================== Tests for single read MutationScorer ============================
//

template <typename R>
class MutationScorerTest : public testing::Test
{
public:
    typedef MutationScorer<R> MS;

protected:
    MutationScorerTest()
        : recursor_(ALL_MOVES, BandingOptions(4, 200)),
          testingConfig_(TestingConfig())
    {}

    virtual ~MutationScorerTest() {}

protected:
    typename MS::RecursorType recursor_;
    QuiverConfig testingConfig_;
};

#define MS MutationScorer<TypeParam>
#define E  typename TypeParam::EvaluatorType

#define recursor  (this->recursor_)
#define params    (this->testingConfig_.QvParams)
#define config    (this->testingConfig_)

extern QvRead AnonymousRead(std::string seq);


TYPED_TEST(MutationScorerTest, BasicTest)
{
    std::string tpl = "GATTACA";
    QvRead read = AnonymousRead("GATTACA");
    E ev(read, tpl, params, true, true);
    MS ms(ev, recursor);
    Mutation mergeableInsertMutation(INSERTION, 4, 'A');
    Mutation unmergeableInsertMutation(INSERTION, 4, 'G');
    Mutation substitutionMutation(SUBSTITUTION, 4, 'T');
    Mutation deletionMutation(DELETION, 4, '-');

    // Testing mutations should not change the template.
    // Let's just make sure of that.
    EXPECT_EQ("GATTACA"        , ms.Template());
    EXPECT_EQ(0                , ms.Score());
    EXPECT_EQ("GATTACA"        , ms.Template());
    EXPECT_EQ(params.Merge[0]  , ms.ScoreMutation(mergeableInsertMutation));
    EXPECT_EQ("GATTACA"        , ms.Template());
    EXPECT_EQ(params.DeletionN , ms.ScoreMutation(unmergeableInsertMutation));
    EXPECT_EQ("GATTACA"        , ms.Template());
    EXPECT_EQ(params.Mismatch  , ms.ScoreMutation(substitutionMutation));
    EXPECT_EQ("GATTACA"        , ms.Template());
    EXPECT_EQ(params.Nce       , ms.ScoreMutation(deletionMutation));
    EXPECT_EQ("GATTACA"        , ms.Template());
}


TYPED_TEST(MutationScorerTest, CopyTest)
{
    std::string tpl = "GATTACA";
    QvRead read = AnonymousRead("GATTACA");
    E ev(read, tpl, params, true, true);
    MS ms(ev, recursor);
    MS msCopy(ms);
    ASSERT_EQ(ms.Score(), msCopy.Score());
}


TYPED_TEST(MutationScorerTest, MutationsAtBeginning)
{
    std::string tpl = "GATTACA";
    QvRead read = AnonymousRead("GATTACA");
    E ev(read, tpl, params, true, true);
    MS ms(ev, recursor);
    Mutation insertBefore(INSERTION, 0, 'A');
    Mutation mergeableInsertMutation1(INSERTION, 1, 'G');
    Mutation mergeableInsertMutation2(INSERTION, 1, 'A');
    Mutation unmergeableInsertMutation(INSERTION, 1, 'T');
    Mutation substitutionMutation(SUBSTITUTION, 0, 'T');
    Mutation deletionMutation(DELETION, 0, '-');

    EXPECT_EQ(0                , ms.Score());
    EXPECT_EQ(params.DeletionN , ms.ScoreMutation(insertBefore));
    EXPECT_EQ(params.Merge[0]  , ms.ScoreMutation(mergeableInsertMutation1));
    EXPECT_EQ(params.Merge[0]  , ms.ScoreMutation(mergeableInsertMutation2));
    EXPECT_EQ(params.DeletionN , ms.ScoreMutation(unmergeableInsertMutation));
    EXPECT_EQ(params.Mismatch  , ms.ScoreMutation(substitutionMutation));
    EXPECT_EQ(params.Nce       , ms.ScoreMutation(deletionMutation));
}

TYPED_TEST(MutationScorerTest, MutationsAtEnd)
{
    std::string tpl = "GATTACA";
    QvRead read = AnonymousRead("GATTACA");
    E ev(read, tpl, params, true, true);
    MS ms(ev, recursor);
    Mutation mergeableInsertMutation(INSERTION, 7, 'A');
    Mutation unmergeableInsertMutation(INSERTION, 7, 'G');
    Mutation substitutionMutation(SUBSTITUTION, 6, 'T');
    Mutation deletionMutation(DELETION, 6, '-');

    EXPECT_EQ(0                , ms.Score());
    EXPECT_EQ(params.Merge[0]  , ms.ScoreMutation(mergeableInsertMutation));
    EXPECT_EQ(params.DeletionN , ms.ScoreMutation(unmergeableInsertMutation));
    EXPECT_EQ(params.Mismatch  , ms.ScoreMutation(substitutionMutation));
    EXPECT_EQ(params.Nce       , ms.ScoreMutation(deletionMutation));
}


TYPED_TEST(MutationScorerTest, TinyTemplate)
{
    std::string tpl = "GTGC";
    QvRead read = AnonymousRead("GTGC");
    E ev(read, tpl, params, true, true);
    MS ms(ev, recursor);

    Mutation deletionAtBeginning(DELETION, 0, '-');
    Mutation deletionAtEnd(DELETION, 3, '-');
    EXPECT_EQ(params.Nce, ms.ScoreMutation(deletionAtBeginning));
    EXPECT_EQ(params.Nce, ms.ScoreMutation(deletionAtEnd));


    Mutation insertAtBeginning(INSERTION, 0, 'T');
    Mutation insertAtEnd(INSERTION, 4, 'T');
    EXPECT_EQ(params.DeletionN, ms.ScoreMutation(insertAtBeginning));
    EXPECT_EQ(params.DeletionN, ms.ScoreMutation(insertAtEnd));

    for (int pos = 0; pos < tpl.length(); pos++)
    {
        Mutation m(SUBSTITUTION, pos, 'A');
        EXPECT_EQ(params.Mismatch, ms.ScoreMutation(m));
    }
}


TYPED_TEST(MutationScorerTest, TemplateMutationWorkflow)
{
    std::string tpl = "GATTACA";
    QvRead read = AnonymousRead("GATTACA");
    E ev(read, tpl, params, true, true);
    MS ms(ev, recursor);
    Mutation insertMutation(INSERTION, 4, 'A');

    EXPECT_EQ("GATTACA"       , ms.Template());
    EXPECT_EQ(0               , ms.Score());
    EXPECT_EQ("GATTACA"       , ms.Template());
    EXPECT_EQ(params.Merge[0] , ms.ScoreMutation(insertMutation));

    std::string newTpl = ApplyMutation(insertMutation, tpl);
    ms.Template(newTpl);
    EXPECT_EQ(params.Merge[0] , ms.Score());
    EXPECT_EQ("GATTAACA"      , ms.Template());
}




TYPED_TEST(MutationScorerTest, DinucleotideInsertionTest)
{
    //                     0123456789012345678
    std::string tplTT   = "CCCCCGATTACACCCCC";
    std::string tplTTTT = "CCCCCGATTTTACACCCCC";
    std::string tplGCTT = "CCCCCGAGCTTACACCCCC";
    std::string tplAATT = "CCCCCGAAATTACACCCCC";

    QvRead read = AnonymousRead("CCCCCGATTTTACACCCCC");
    ReadScorer ez(config);
    float scoreTTTT = ez.Score(tplTTTT, read);
    EXPECT_EQ(0, scoreTTTT);

    QvEvaluator e(read, tplTT, params);
    SparseSimpleQvRecursor r(config.MovesAvailable, config.Banding);
    SparseSimpleQvMutationScorer ms(e, r);

    EXPECT_EQ(0, ms.ScoreMutation(Mutation(INSERTION, 7, 7, "TT")));
    EXPECT_EQ(0, ms.ScoreMutation(Mutation(INSERTION, 8, 8, "TT")));
    EXPECT_EQ(0, ms.ScoreMutation(Mutation(INSERTION, 9, 9, "TT")));

    EXPECT_EQ(ez.Score(tplGCTT, read), ms.ScoreMutation(Mutation(INSERTION, 7, 7, "GC")));
    EXPECT_EQ(ez.Score(tplAATT, read), ms.ScoreMutation(Mutation(INSERTION, 7, 7, "AA")));
    EXPECT_EQ(ez.Score(tplAATT, read), ms.ScoreMutation(Mutation(INSERTION, 6, 6, "AA")));
}

TYPED_TEST(MutationScorerTest, DinucleotideDeletionTest)
{
    //                     0123456789012345678
    std::string tplTT   = "CCCCCGATTACACCCCC";
    std::string tplTTTT = "CCCCCGATTTTACACCCCC";
    std::string tplGCTT = "CCCCCGAGCTTACACCCCC";

    QvRead read = AnonymousRead("CCCCCGATTACACCCCC");
    ReadScorer ez(config);
    float scoreTT = ez.Score(tplTT, read);
    EXPECT_EQ(0, scoreTT);

    QvEvaluator e(read, tplTTTT, params);
    SparseSimpleQvRecursor r(config.MovesAvailable, config.Banding);
    SparseSimpleQvMutationScorer ms(e, r);

    EXPECT_EQ(scoreTT, ms.ScoreMutation(Mutation(DELETION, 7, 9, "")));
    EXPECT_EQ(scoreTT, ms.ScoreMutation(Mutation(DELETION, 8, 10, "")));
    EXPECT_EQ(scoreTT, ms.ScoreMutation(Mutation(DELETION, 9, 11, "")));

    QvEvaluator e2(read, tplGCTT, params);
    SparseSimpleQvRecursor r2(config.MovesAvailable, config.Banding);
    SparseSimpleQvMutationScorer ms2(e2, r2);
    EXPECT_EQ(scoreTT, ms2.ScoreMutation(Mutation(DELETION, 7, 9, "")));
}
