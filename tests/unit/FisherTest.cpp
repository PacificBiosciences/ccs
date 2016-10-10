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

// Author: Lance Hepler

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <pacbio/statistics/Fisher.h>

using std::string;

using namespace PacBio::Statistics;  // NOLINT

namespace {

TEST(FisherTest, EquivalentWithR)
{
    // sapply(0:100, function(x) { fisher.test(matrix(c(x,1000,10,1000), nrow=2, byrow = TRUE), alternative = "greater")$p.value})
    std::vector<double> pValuesFromR{
        1.000000e+00, 9.995009e-01, 9.967696e-01, 9.886124e-01, 9.710086e-01, 9.403090e-01,
        8.943966e-01, 8.333319e-01, 7.593288e-01, 6.761990e-01, 5.885378e-01, 5.009200e-01,
        4.172848e-01, 3.405813e-01, 2.726612e-01, 2.143576e-01, 1.656755e-01, 1.260236e-01,
        9.444223e-02, 6.979418e-02, 5.090970e-02, 3.668358e-02, 2.613167e-02, 1.841604e-02,
        1.284825e-02, 8.879206e-03, 6.081771e-03, 4.130825e-03, 2.783570e-03, 1.861740e-03,
        1.236419e-03, 8.156574e-04, 5.346868e-04, 3.484045e-04, 2.257325e-04, 1.454637e-04,
        9.325708e-05, 5.949549e-05, 3.778008e-05, 2.388437e-05, 1.503578e-05, 9.427211e-06,
        5.887959e-06, 3.663910e-06, 2.271922e-06, 1.404031e-06, 8.648820e-07, 5.311204e-07,
        3.251928e-07, 1.985425e-07, 1.208873e-07, 7.341252e-08, 4.446997e-08, 2.687285e-08,
        1.620137e-08, 9.745838e-09, 5.849967e-09, 3.504205e-09, 2.094890e-09, 1.249972e-09,
        7.444536e-10, 4.425907e-10, 2.626772e-10, 1.556412e-10, 9.207351e-11, 5.438474e-11,
        3.207564e-11, 1.889088e-11, 1.111038e-11, 6.525696e-12, 3.827944e-12, 2.242665e-12,
        1.312326e-12, 7.670354e-13, 4.478190e-13, 2.611677e-13, 1.521536e-13, 8.855353e-14,
        5.148791e-14, 2.990849e-14, 1.735751e-14, 1.006462e-14, 5.830912e-15, 3.375342e-15,
        1.952329e-15, 1.128375e-15, 6.516743e-16, 3.760927e-16, 2.168987e-16, 1.250049e-16,
        7.199700e-17, 4.144087e-17, 2.383854e-17, 1.370489e-17, 7.874537e-18, 4.522064e-18,
        2.595485e-18, 1.488944e-18, 8.537362e-19, 4.892853e-19, 2.802857e-19};
    for (int i = 11; i < 100; ++i)
        EXPECT_NEAR(pValuesFromR.at(i), Fisher::fisher_exact_tiss(i, 1000, 10, 1000),
                    pValuesFromR.at(i) / 1e5);
}
}
