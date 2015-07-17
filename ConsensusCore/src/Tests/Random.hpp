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

#pragma once

#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <iostream>
#include <string>
#include <vector>

#include <ConsensusCore/Mutation.hpp>
#include <ConsensusCore/Features.hpp>
#include <ConsensusCore/Quiver/QvEvaluator.hpp>
#include <ConsensusCore/Read.hpp>
#include <ConsensusCore/Utils.hpp>

#include "ParameterSettings.hpp"

using ConsensusCore::Mutation;
using ConsensusCore::MutationType;
using ConsensusCore::QvEvaluator;
using ConsensusCore::QvModelParams;
using ConsensusCore::QvSequenceFeatures;
using ConsensusCore::QvRead;

template<typename RNG>
std::string
RandomSequence(RNG& rng, int length)
{
    const char* bases = "ACGT";
    boost::random::uniform_int_distribution<> indexDist(0, 3);
    std::stringstream ss;
    for (int i = 0; i < length; ++i)
    {
        ss << bases[indexDist(rng)];
    }
    return ss.str();
}

template<typename RNG>
float*
RandomQvArray(RNG& rng, int length)
{
    float* array = new float[length];
    boost::random::uniform_int_distribution<> indexDist(0, 127);
    for (int i = 0; i < length; ++i)
    {
        array[i] = static_cast<float>(indexDist(rng));
    }
    return array;
}

template<typename RNG>
float*
RandomTagArray(RNG& rng, int length)
{
    std::string seq = RandomSequence(rng, length);
    float* array = new float[length];
    for (int i = 0; i < length; ++i)
    {
        array[i] = static_cast<float>(seq[i]);
    }
    return array;
}

template<typename RNG>
int
RandomPoissonDraw(RNG& rng, int mean)
{
    boost::random::poisson_distribution<> dist(mean);
    return dist(rng);
}

template<typename RNG>
bool
RandomBernoulliDraw(RNG& rng, float p)
{
    boost::random::bernoulli_distribution<> dist(p);
    return dist(rng);
}

template<typename RNG>
QvEvaluator
RandomQvEvaluator(RNG& rng, int length)
{
    std::string tpl = RandomSequence(rng, length);

    int readLength = RandomPoissonDraw(rng, length);
    std::string seq = RandomSequence(rng, readLength);

    float* insQv = RandomQvArray(rng, readLength);
    float* subsQv = RandomQvArray(rng, readLength);
    float* delQv = RandomQvArray(rng, readLength);
    float* delTag = RandomTagArray(rng, readLength);
    float* mergeQv = RandomQvArray(rng, readLength);

    QvSequenceFeatures f(seq, insQv, subsQv, delQv, delTag, mergeQv);
    QvRead read(f, "anonymous", "unknown");

    delete[] insQv;
    delete[] subsQv;
    delete[] delQv;
    delete[] delTag;
    delete[] mergeQv;

    bool pinStart = RandomBernoulliDraw(rng, 0.5);
    bool pinEnd = RandomBernoulliDraw(rng, 0.5);
    return QvEvaluator(read, tpl, TestingParams(), pinStart, pinEnd);
}

template<typename RNG>
std::vector<int>
RandomSampleWithoutReplacement(RNG& rng, int n, int k)
{
    // Random sample of k elements from [0..n) without replacement
    std::vector<int> draws;
    boost::random::uniform_int_distribution<> indexDist(0, n - 1);
    while (draws.size() < k) {
        int draw = indexDist(rng);
        if (std::find(draws.begin(), draws.end(), draw) == draws.end())
        {
            draws.push_back(draw);
        }
    }
}

typedef boost::mt19937 Rng;
