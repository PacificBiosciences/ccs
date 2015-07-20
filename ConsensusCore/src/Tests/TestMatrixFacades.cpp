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

#include <algorithm>
#include <boost/tuple/tuple.hpp>
#include <climits>
#include <iostream>
#include <string>
#include <typeinfo>
#include <vector>

#include <ConsensusCore/LValue.hpp>
#include <ConsensusCore/Matrix/DenseMatrix.hpp>
#include <ConsensusCore/Matrix/SparseMatrix.hpp>

using std::cout;
using std::endl;

using ConsensusCore::DenseMatrixF;
using ConsensusCore::SparseMatrixF;
using ConsensusCore::lfloat;


template <typename T>
class MatrixTest : public ::testing::Test
{
public:
    virtual ~MatrixTest() {}
};

using testing::Types;
// typedef Types<DenseMatrix> Implementations;
typedef Types<DenseMatrixF, SparseMatrixF> Implementations;
TYPED_TEST_CASE(MatrixTest, Implementations);


TYPED_TEST(MatrixTest, Basic)
{
    float NEG_INF = -FLT_MAX;
    EXPECT_EQ(NEG_INF, lfloat());

    TypeParam m(10, 10);
    EXPECT_EQ(10, m.Rows());
    EXPECT_EQ(10, m.Columns());
    for (int i = 0; i < 10; ++i)
    {
        for (int j = 0; j < 10; ++j)
        {
            EXPECT_EQ(NEG_INF, m(i, j));
        }
    }

    m.StartEditingColumn(1, 0, 10);
    m.Set(1, 1, 5);
    m.Set(2, 1, 6);
    m.FinishEditingColumn(1, 1, 3);
    EXPECT_EQ(5, m(1, 1));
    EXPECT_EQ(6, m(2, 1));
    m.ClearColumn(1);
    for (int j = 0; j < 10; ++j)
    {
        EXPECT_EQ(NEG_INF, m(j, 1));
    }
}

TYPED_TEST(MatrixTest, Nullability)
{
    TypeParam m(10, 10);
    EXPECT_TRUE(!m.IsNull());
    EXPECT_TRUE(TypeParam::Null().IsNull());
    TypeParam nullCopy = TypeParam::Null();
    EXPECT_TRUE(nullCopy.IsNull());
}

TYPED_TEST(MatrixTest, Ranges)
{
    TypeParam m(10, 10);

    // check settings on empty matrix
    EXPECT_EQ(0, m.UsedEntries());
    for (int j = 0; j < 10; j++)
    {
        int start, end;
        boost::tie(start, end) = m.UsedRowRange(j);
        EXPECT_EQ(0, start);
        EXPECT_EQ(0, end);
    }

    // simple modifications
    for (int j = 0; j < 10; j++)
    {
        int start, end;
        m.StartEditingColumn(j, 0, 10);
        m.Set(2, j, 0.0);
        m.Set(3, j, 0.0);
        m.Set(4, j, 0.0);
        m.FinishEditingColumn(j, 2, 5);
        boost::tie(start, end) = m.UsedRowRange(j);
        EXPECT_EQ(2, start);
        EXPECT_EQ(5, end);
    }
    EXPECT_EQ(30, m.UsedEntries());

    for (int j = 0; j < 10; j++)
    {
        int start, end;
        m.ClearColumn(j);
        boost::tie(start, end) = m.UsedRowRange(j);
        EXPECT_EQ(0, start);
        EXPECT_EQ(0, end);
    }
    EXPECT_EQ(0, m.UsedEntries());
}

TYPED_TEST(MatrixTest, IsColumnEmpty)
{
    TypeParam m(10, 10);

    EXPECT_EQ(m.IsColumnEmpty(0), true);

    m.StartEditingColumn(0, 0, 0);
    m.Set(1, 0, 0.0);
    m.FinishEditingColumn(0, 0, 5);

    EXPECT_EQ(m.IsColumnEmpty(0), false);
}

TYPED_TEST(MatrixTest, SSE)
{
    TypeParam m(10, 10);
    const float cookieArray[] = {0, 1, 2, 3};
    float cookieReadArray[4];
    __m128 cookie = _mm_loadu_ps(cookieArray);

    // test Set4
    m.StartEditingColumn(0, 0, 0);
    m.Set4(0, 0, cookie);
    m.FinishEditingColumn(0, 0, 4);
    EXPECT_EQ(0, m(0, 0));
    EXPECT_EQ(1, m(1, 0));
    EXPECT_EQ(2, m(2, 0));
    EXPECT_EQ(3, m(3, 0));

    // test Get4
    __m128 cookieRead = m.Get4(0, 0);
    _mm_storeu_ps(cookieReadArray, cookieRead);
    for (int i = 0; i < 4; ++i)
    {
        EXPECT_EQ(cookieArray[i], cookieReadArray[i]);
    }
}

TYPED_TEST(MatrixTest, ToHostArray)
{
    TypeParam m(10, 10);
    int v = 0;
    for (int j = 0; j < 10; j++)
    {
        m.StartEditingColumn(j, 0, 0);
        for (int i = 0; i < 10; i++)
        {
            m.Set(i, j, v);
            v++;
        }
        m.FinishEditingColumn(j, 0, 10);
    }

    float* hostArray;
    int rows, cols;
    m.ToHostMatrix(&hostArray, &rows, &cols);
    EXPECT_EQ(10, rows);
    EXPECT_EQ(10, cols);
    v = 0;
    for (int j = 0; j < 10; j++)
    {
        for (int i = 0; i < 10; i++)
        {
            EXPECT_EQ(v, hostArray[i * cols + j]);
            v++;
        }
    }
    delete[] hostArray;
}


TYPED_TEST(MatrixTest, NonSequentialAccess)
{}


TYPED_TEST(MatrixTest, Holes)
{
    // Test inserting -FLT_MAX values into the middle of the matrix,
    // make sure it doesn't confuse the matrix class
}

TYPED_TEST(MatrixTest, BigBandedMatrix)
{
    // fill a big matrix with constant bandwidth
    const int bandWidth = 5;
    const int M = 1000;
    const int N = 1000;
    TypeParam m(M, N);
    for (int j = 0; j < N; j++)
    {
        m.StartEditingColumn(j, 0, 0);
        int start = std::max(0, j - bandWidth);
        int end   = std::min(M, j + bandWidth + 1);
        for (int i = start; i < end; i++)
        {
            m.Set(i, j, i / (1. + j));
        }
        m.FinishEditingColumn(j, start, end);
    }

    cout << typeid(m).name() << " : " << m.AllocatedEntries() << endl;
}

TYPED_TEST(MatrixTest, BigIrregularBandedMatrix)
{
    // fill a big matrix, experimenting with modulating the bandwidth
}


TYPED_TEST(MatrixTest, CopyTest)
{
    TypeParam m(4, 4);
    m.StartEditingColumn(1, 0, 4);
    m.Set(1, 1, 5);
    m.FinishEditingColumn(1, 1, 2);
    TypeParam mCopy(m);

    ASSERT_EQ(5, mCopy(1, 1));
}
