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

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <iostream>
#include <string>
#include <vector>

#include <ConsensusCore/Quiver/QvEvaluator.hpp>
#include <ConsensusCore/Features.hpp>
#include <ConsensusCore/Utils.hpp>

#include "Random.hpp"
#include "SseTestingUtils.hpp"
#include "ParameterSettings.hpp"

using namespace ConsensusCore; // NOLINT
using std::cout;
using std::endl;

class QvEvaluatorTest : public testing::Test
{
public:
    void SetUp()
    {
        int numEvaluators = 200;
        int tplLen = 20;
        int seed = 42;

        Rng rng(seed);
        for (int n = 0; n < numEvaluators; n++)
        {
            fuzzEvaluators_.push_back(RandomQvEvaluator(rng, tplLen));
        }
    }

    virtual ~QvEvaluatorTest() {}

protected:
    std::vector<QvEvaluator> fuzzEvaluators_;
};


TEST_F(QvEvaluatorTest, IncVsInc4)
{
    foreach (const QvEvaluator& e, this->fuzzEvaluators_)
    {
        int I = e.ReadLength();
        int J = e.TemplateLength();

        for (int j = 0; j <= J - 1 ; j++)
            for (int i = 0; i <= I - 4; i++)
                COMPARE4(e.Inc4, e.Inc, i, j);
    }
}


TEST_F(QvEvaluatorTest, DelVsDel4)
{
    foreach (const QvEvaluator& e, this->fuzzEvaluators_)
    {
        int I = e.ReadLength();
        int J = e.TemplateLength();

        for (int j = 0; j <= J - 1 ; j++)
            for (int i = 0; i <= I - 3; i++)
                COMPARE4(e.Del4, e.Del, i, j);
    }
}


TEST_F(QvEvaluatorTest, ExtraVsExtra4)
{
    foreach (const QvEvaluator& e, this->fuzzEvaluators_)
    {
        int I = e.ReadLength();
        int J = e.TemplateLength();

        for (int j = 0; j <= J ; j++)
            for (int i = 0; i <= I - 4; i++)
                COMPARE4(e.Extra4, e.Extra, i, j);
    }
}


TEST_F(QvEvaluatorTest, MergeVsMerge4)
{
    foreach (const QvEvaluator& e, this->fuzzEvaluators_)
    {
        int I = e.ReadLength();
        int J = e.TemplateLength();

        for (int j = 0; j <= J - 2; j++)
            for (int i = 0; i <= I - 4; i++)
                COMPARE4(e.Merge4, e.Merge, i, j);
    }
}


TEST_F(QvEvaluatorTest, BadTagTest)
{
    Rng rng(42);
    int n = 20;

    int readLength = RandomPoissonDraw(rng, n);
    std::string seq = RandomSequence(rng, readLength);

    float* insQv = RandomQvArray(rng, readLength);
    float* subsQv = RandomQvArray(rng, readLength);
    float* delQv = RandomQvArray(rng, readLength);
    float* delTag = RandomTagArray(rng, readLength);
    float* mergeQv = RandomQvArray(rng, readLength);
    delTag[5] = 3;

    QvSequenceFeatures* f;
    ASSERT_THROW(f = new QvSequenceFeatures(seq, insQv, subsQv, delQv, delTag, mergeQv),
                 InternalError);

    delete[] insQv;
    delete[] subsQv;
    delete[] delQv;
    delete[] delTag;
    delete[] mergeQv;
}
