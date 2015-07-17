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

#include <boost/format.hpp>
#include <iostream>
#include <string>
#include <vector>

#include <ConsensusCore/Align/PairwiseAlignment.hpp>
#include <ConsensusCore/Matrix/DenseMatrix.hpp>
#include <ConsensusCore/Matrix/SparseMatrix.hpp>
#include <ConsensusCore/Quiver/QvEvaluator.hpp>
#include <ConsensusCore/Quiver/QuiverConfig.hpp>
#include <ConsensusCore/Quiver/SimpleRecursor.hpp>
#include <ConsensusCore/Quiver/SseRecursor.hpp>
#include <ConsensusCore/Features.hpp>

#include "MatrixPrinting.hpp"
#include "ParameterSettings.hpp"
#include "Random.hpp"

using namespace ConsensusCore; // NOLINT

//
//  Instantiate the concrete test classes, by speciying the implementations we seek to test.
//
typedef testing::Types<SimpleQvRecursor,
                       SseQvRecursor,
                       SparseSimpleQvRecursor,
                       SparseSseQvRecursor> Implementations;

TYPED_TEST_CASE(RecursorTest    , Implementations);
TYPED_TEST_CASE(RecursorFuzzTest, Implementations);


template <typename T>
class RecursorTest : public testing::Test
{
protected:
    RecursorTest()
        : noBanding_(0, 1e9),
          standardBanding_(4, 200),
          testingParams_(TestingParams())
    {}

    virtual ~RecursorTest() {}

protected:
    BandingOptions noBanding_;
    BandingOptions standardBanding_;
    typename T::EvaluatorType::ParamsType testingParams_;
};

// Macros to simplify type lookup code.
#define R TypeParam
#define M typename TypeParam::MatrixType
#define E typename TypeParam::EvaluatorType
#define NULL_MATRIX TypeParam::MatrixType::Null()


extern QvRead AnonymousRead(std::string seq);

TYPED_TEST(RecursorTest, SmallMatchTest)
{
    std::string tpl("GATG");
    QvRead read = AnonymousRead("GATG");
    E e(read, tpl, this->testingParams_);

    R recursor(BASIC_MOVES, this->noBanding_);
    M alpha(tpl.length() + 1, read.Length() + 1);
    M beta(tpl.length() + 1, read.Length() + 1);
    recursor.FillAlpha(e, NULL_MATRIX, alpha);
    recursor.FillBeta(e, NULL_MATRIX, beta);

    EXPECT_FLOAT_EQ(0.0f, alpha(read.Length(), tpl.length()));
    EXPECT_FLOAT_EQ(0.0f, beta(0, 0));
    //    std::cout << std::endl;
    //    std::cout << ConsensusCore::detail::PrintMatrix(alpha) << std::endl;

    // try the traceback
    const PairwiseAlignment* alignment = recursor.Alignment(e, alpha);
    EXPECT_EQ("GATG", alignment->Target());
    EXPECT_EQ("GATG", alignment->Query());
    delete alignment;

    // Make sure Beta gave the same score
    EXPECT_FLOAT_EQ(0.0f, beta(0, 0));
    //    std::cout << std::endl;
    //    std::cout << ConsensusCore::detail::PrintMatrix(beta) << std::endl;
}

TYPED_TEST(RecursorTest, SmallMismatchTest)
{
    std::string tpl("GATG");
    QvRead read = AnonymousRead("GATC");
    E e(read, tpl, this->testingParams_);

    R recursor(BASIC_MOVES, this->noBanding_);
    M alpha(read.Length() + 1, tpl.length() + 1);
    M beta(read.Length() + 1, tpl.length() + 1);
    recursor.FillAlpha(e, NULL_MATRIX, alpha);
    recursor.FillBeta(e, NULL_MATRIX, beta);

    EXPECT_FLOAT_EQ(-10.0f, alpha(read.Length(), tpl.length()))
        << std::endl << PrintMatrix(alpha) << std::endl;

    EXPECT_FLOAT_EQ(-10.0f, beta(0, 0))
        << std::endl << PrintMatrix(alpha) << std::endl;

    // try the traceback
    const PairwiseAlignment* alignment = recursor.Alignment(e, alpha);
    EXPECT_EQ("GATG", alignment->Target());
    EXPECT_EQ("GATC", alignment->Query());
    delete alignment;
}


TYPED_TEST(RecursorTest, SmallMergeTest)
{
    std::string tpl("GATT");
    QvRead read = AnonymousRead("GAT");
    E e(read, tpl, this->testingParams_);

    R recursor(BASIC_MOVES | MERGE, this->noBanding_);
    M alpha(read.Length() + 1, tpl.length() + 1);
    M beta(read.Length() + 1, tpl.length() + 1);
    recursor.FillAlpha(e, NULL_MATRIX, alpha);
    recursor.FillBeta(e, NULL_MATRIX, beta);

    EXPECT_FLOAT_EQ(-2.0f, alpha(read.Length(), tpl.length()));
    EXPECT_FLOAT_EQ(-2.0f, beta(0, 0));
    //    std::cout << std::endl;
    //    std::cout << ConsensusCore::detail::PrintMatrix(alpha) << std::endl;

    // try the traceback
    const PairwiseAlignment* alignment = recursor.Alignment(e, alpha);
    EXPECT_EQ("GATT", alignment->Target());
    EXPECT_EQ("GA-T", alignment->Query());
    delete alignment;

    // Make sure Beta gave the same score
    EXPECT_FLOAT_EQ(-2.0f, beta(0, 0));
    //    std::cout << std::endl;
    //    std::cout << ConsensusCore::detail::PrintMatrix(beta) << std::endl;
}


TYPED_TEST(RecursorTest, MediumSizedTest)
{
    // This is the same as the test in demo.py.
    //    tpl  = "GATTACA"*10
    //    read = "GATTACA"*3 + "GATTTTTTACA"*4 + "GATTACA"*3
    std::stringstream tplSs, readSs;
    for (int i = 0; i < 10; i++) tplSs  << "GATTACA";
    for (int i = 0; i < 3; i++)  readSs << "GATTACA";
    for (int i = 3; i < 7; i++)  readSs << "GATTTTTTACA";
    for (int i = 7; i < 10; i++) readSs << "GATTACA";

    std::string tpl(tplSs.str());
    QvRead read = AnonymousRead(readSs.str());
    E e(read, tpl, this->testingParams_);

    R recursor(BASIC_MOVES, this->standardBanding_);
    M alpha(read.Length() + 1, tpl.length() + 1);
    M beta(read.Length() + 1, tpl.length() + 1);
    recursor.FillAlphaBeta(e, alpha, beta);
    EXPECT_FLOAT_EQ(-80.0f, beta(0, 0));
}

TYPED_TEST(RecursorTest, LinkTest)
{
    std::string tpl("GATTCTC");
    QvRead read = AnonymousRead("GATCTTC");
    E e(read, tpl, this->testingParams_);

    R recursor(BASIC_MOVES | MERGE, this->noBanding_);
    M alpha(read.Length() + 1, tpl.length() + 1);
    M beta(read.Length() + 1, tpl.length() + 1);
    recursor.FillAlphaBeta(e, alpha, beta);

    //    std::cout << std::endl;
    //    std::cout << ConsensusCore::detail::PrintMatrix(alpha) << std::endl;
    //    std::cout << std::endl;
    //    std::cout << ConsensusCore::detail::PrintMatrix(beta) << std::endl;

    float score = beta(0, 0);
    for (int j = 2; j < tpl.length() - 1; j++)
    {
        float linkScore = recursor.LinkAlphaBeta(e, alpha, j, beta, j, j);
        ASSERT_FLOAT_EQ(score, linkScore)
            << "(Column " << j << ")";
    }
}


TYPED_TEST(RecursorTest, ExtendAlphaTest)
{
    std::string tpl("GATTCTC");
    QvRead read = AnonymousRead("GATCTTC");
    E e(read, tpl, this->testingParams_);

    R recursor(BASIC_MOVES | MERGE, this->noBanding_);
    M alpha(read.Length() + 1, tpl.length() + 1);
    M beta(read.Length() + 1, tpl.length() + 1);
    recursor.FillAlphaBeta(e, alpha, beta);

    //    std::cout << std::endl;
    //    std::cout << ConsensusCore::detail::PrintMatrix(alpha) << std::endl;
    //    std::cout << std::endl;
    //    std::cout << ConsensusCore::detail::PrintMatrix(beta) << std::endl;

    M ext(read.Length() + 1, 2);
    for (int j = 2; j <= tpl.length() - 1; j++)
    {
        recursor.ExtendAlpha(e, alpha, j, ext);
        for (int extCol = 0; extCol < 2; extCol++)
        {
            for (int i = 0; i <= read.Length(); i++)
            {
                ASSERT_FLOAT_EQ(alpha(i, j + extCol), ext(i, extCol))
                        << i << " " << j << " " << extCol << std::endl;
            }
        }
    }
}


TYPED_TEST(RecursorTest, ExtendBetaTest)
{
    std::string tpl("GATTCTC");
    QvRead read = AnonymousRead("GATCTTC");
    E e(read, tpl, this->testingParams_);

    R recursor(BASIC_MOVES | MERGE, this->noBanding_);
    M alpha(read.Length() + 1, tpl.length() + 1);
    M beta(read.Length() + 1, tpl.length() + 1);
    recursor.FillAlphaBeta(e, alpha, beta);

    M ext(read.Length() + 1, 2);
    for (int j = 1; j <= tpl.length() - 2; j++)
    {
        recursor.ExtendBeta(e, beta, j, ext);
        for (int i = 0; i <= read.Length(); i++)
        {
            ASSERT_FLOAT_EQ(beta(i, j), ext(i, 1))
                << i << " " << j << std::endl;
            ASSERT_FLOAT_EQ(beta(i, j - 1), ext(i, 0))
                << i << " " << j-1 << std::endl;
        }
    }
}


// ----------------------------------------------------------------------------
// Fuzz tests --- testing applied to several hundred random templates and reads
// with a goal of catching rare bugs not captured by existing test cases.  Not
// a substitute for focused regression tests!
// ----------------------------------------------------------------------------


template <typename T>
class RecursorFuzzTest : public testing::Test
{
protected:
    RecursorFuzzTest()
        : banding_(4, 200),
          testingParams_(TestingParams())
    {}

    void SetUp()
    {
        int numEvaluators = 200;
        int tplLen = 20;

        Rng rng(42);
        for (int n = 0; n < numEvaluators; n++)
        {
            fuzzEvaluators_.push_back(RandomQvEvaluator(rng, tplLen));
        }
    }

    virtual ~RecursorFuzzTest() {}

protected:
    BandingOptions banding_;
    typename T::EvaluatorType::ParamsType testingParams_;
    std::vector<QvEvaluator> fuzzEvaluators_;
};


TYPED_TEST(RecursorFuzzTest, AlphaBetaConcordance)
{
    R recursor(BASIC_MOVES | MERGE, this->banding_);

    foreach (const QvEvaluator& e, this->fuzzEvaluators_)
    {
        int tplLength = e.TemplateLength();
        int readLength = e.ReadLength();

        M alpha(readLength + 1, tplLength + 1);
        M beta(readLength + 1, tplLength + 1);

        recursor.FillAlphaBeta(e, alpha, beta);
        EXPECT_FLOAT_EQ(alpha(readLength, tplLength), beta(0, 0));
    }
}


TYPED_TEST(RecursorFuzzTest, Alignment)
{
    R recursor(BASIC_MOVES | MERGE, this->banding_);

    foreach (const QvEvaluator& e, this->fuzzEvaluators_)
    {
        int tplLength = e.TemplateLength();
        int readLength = e.ReadLength();

        M alpha(readLength + 1, tplLength + 1);
        M beta(readLength + 1, tplLength + 1);

        recursor.FillAlphaBeta(e, alpha, beta);
        const PairwiseAlignment* alignment = recursor.Alignment(e, alpha);
        EXPECT_TRUE(alignment->Target().length() == alignment->Query().length());
        delete alignment;
    }
}

TYPED_TEST(RecursorFuzzTest, LinkAlphaBeta)
{
    R recursor(BASIC_MOVES | MERGE, this->banding_);

    foreach (const QvEvaluator& e, this->fuzzEvaluators_)
    {
        int tplLength = e.TemplateLength();
        int readLength = e.ReadLength();

        M alpha(readLength + 1, tplLength + 1);
        M beta(readLength + 1, tplLength + 1);

        recursor.FillAlphaBeta(e, alpha, beta);
        float score = beta(0, 0);
        for (int j = 2; j < tplLength - 1; j++)
        {
            float linkScore = recursor.LinkAlphaBeta(e, alpha, j, beta, j, j);
            ASSERT_FLOAT_EQ(score, linkScore)
                << "(Column " << j << ")";
        }
    }
}

TYPED_TEST(RecursorFuzzTest, ExtendAlpha)
{
    R recursor(BASIC_MOVES | MERGE, this->banding_);

    foreach (const QvEvaluator& e, this->fuzzEvaluators_)
    {
        int tplLength = e.TemplateLength();
        int readLength = e.ReadLength();

        M alpha(readLength + 1, tplLength + 1);
        M beta(readLength + 1, tplLength + 1);
        M ext(readLength + 1, 2);

        for (int j = 2; j <= tplLength - 1; j++)
        {
            recursor.ExtendAlpha(e, alpha, j, ext);
            for (int extCol = 0; extCol < 2; extCol++)
            {
                for (int i = 0; i <= readLength; i++)
                {
                    ASSERT_FLOAT_EQ(alpha(i, j + extCol), ext(i, extCol))
                            << i << " " << j << " " << extCol << std::endl;
                }
            }
        }
    }
}

TYPED_TEST(RecursorFuzzTest, ExtendBeta)
{
    R recursor(BASIC_MOVES | MERGE, this->banding_);

    foreach (const QvEvaluator& e, this->fuzzEvaluators_)
    {
        int tplLength = e.TemplateLength();
        int readLength = e.ReadLength();

        M alpha(readLength + 1, tplLength + 1);
        M beta(readLength + 1, tplLength + 1);
        M ext(readLength + 1, 2);

        for (int j = 1; j <= tplLength - 2; j++)
        {
            recursor.ExtendBeta(e, beta, j, ext);
            for (int extCol = 0; extCol < 2; extCol++)
            {
                int jj = j + extCol - 1;
                for (int i = 0; i <= readLength; i++)
                {
                    ASSERT_FLOAT_EQ(beta(i, jj), ext(i, extCol))
                        << i << " " << jj << " " << extCol << std::endl;
                }
            }
        }
    }
}
