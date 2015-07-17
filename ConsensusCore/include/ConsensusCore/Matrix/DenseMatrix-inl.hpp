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

#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/tuple/tuple.hpp>
#include <cassert>

#include <ConsensusCore/Interval.hpp>
#include <ConsensusCore/Utils.hpp>


using boost::numeric::ublas::matrix_column;

namespace ConsensusCore {
    //
    // Nullability
    //
    template<typename F, typename Z>
    inline const DenseMatrix<F, Z>&
    DenseMatrix<F, Z>::Null()
    {
        static DenseMatrix<F, Z>* nullObj = new DenseMatrix<F, Z>(0, 0);
        return *nullObj;
    }

    template<typename F, typename Z>
    inline bool
    DenseMatrix<F, Z>::IsNull() const
    {
        return (Rows() == 0 && Columns() == 0);
    }

    //
    // Size information
    //
    template<typename F, typename Z>
    inline const int
    DenseMatrix<F, Z>::Rows() const
    {
        return boost_dense_matrix::size1();
    }

    template<typename F, typename Z>
    inline const  int
    DenseMatrix<F, Z>::Columns() const
    {
        return boost_dense_matrix::size2();
    }

    //
    // Entry range queries per column
    //
    template<typename F, typename Z>
    inline void
    DenseMatrix<F, Z>::StartEditingColumn(int j, int hintBegin, int hintEnd)
    {
        assert(columnBeingEdited_ == -1);
        columnBeingEdited_ = j;
        ClearColumn(j);
    }

    template<typename F, typename Z>
    inline void
    DenseMatrix<F, Z>::FinishEditingColumn(int j, int usedRowsBegin, int usedRowsEnd)
    {
        assert(columnBeingEdited_ == j);
        usedRanges_[j] = Interval(usedRowsBegin, usedRowsEnd);
        DEBUG_ONLY(CheckInvariants(columnBeingEdited_));
        columnBeingEdited_ = -1;
    }

    template<typename F, typename Z>
    inline Interval
    DenseMatrix<F, Z>::UsedRowRange(int j) const
    {
        assert(0 <= j && j < (int)usedRanges_.size());
        return usedRanges_[j];
    }

    template<typename F, typename Z>
    inline bool
    DenseMatrix<F, Z>::IsColumnEmpty(int j) const
    {
        assert(0 <= j && j < (int)usedRanges_.size());
        return (usedRanges_[j].Begin >= usedRanges_[j].End);
    }

    //
    // Accessors
    //
    template<typename F, typename Z>
    inline void
    DenseMatrix<F, Z>::Set(int i, int j, F v)
    {
        assert(columnBeingEdited_ == j);
        boost_dense_matrix::operator()(i, j) = v;
    }

    template<typename F, typename Z>
    inline bool
    DenseMatrix<F, Z>::IsAllocated(int i, int j) const
    {
        assert(0 <= i && i < Rows() && 0 <= j && j < Columns());
        return true;
    }

    template<typename F, typename Z>
    inline F
    DenseMatrix<F, Z>::Get(int i, int j) const
    {
        return (*this)(i, j);
    }

    template<typename F, typename Z>
    inline const F&
    DenseMatrix<F, Z>::operator() (int i, int j) const
    {
        return boost_dense_matrix::operator()(i, j);
    }

    template<typename F, typename Z>
    inline void
    DenseMatrix<F, Z>::ClearColumn(int j)
    {
        DEBUG_ONLY(CheckInvariants(j);)
        // (Rely on the fact that the underlying memory is stored
        // contiguously)
        int begin, end;
        boost::tie(begin, end) = usedRanges_[j];
        std::fill_n((F*)&boost_dense_matrix::operator()(begin, j),  // NOLINT
                    end - begin,
                    typename boost_dense_matrix::value_type());
        usedRanges_[j] = Interval(0, 0);
        DEBUG_ONLY(CheckInvariants(j);)
    }

    //
    // SSE
    //
    template<>
    inline __m128
    DenseMatrix<float, lvalue<float>>::Get4(int i, int j) const
    {
        assert(0 <= i && i <= Rows() - 4);
        return _mm_loadu_ps(&boost_dense_matrix::operator()(i, j).value);
    }

    template<typename F, typename Z>
    inline __m128
    DenseMatrix<F, Z>::Get4(int i, int j) const
    {
        throw std::runtime_error("cannot use Get4 with non-f32 type!");
    }

    template<>
    inline void
    DenseMatrix<float, lvalue<float>>::Set4(int i, int j, __m128 v4)
    {
        assert(columnBeingEdited_ == j);
        assert(0 <= i && i <= Rows() - 4);
        _mm_storeu_ps(&boost_dense_matrix::operator()(i, j).value, v4);
    }

    template<typename F, typename Z>
    inline void
    DenseMatrix<F, Z>::Set4(int i, int j, __m128 v4)
    {
        throw std::runtime_error("cannot use Set4 with non-f32 type!");
    }

    //
    // Performance insensitive routines are not inlined
    //
    template<typename F, typename Z>
    DenseMatrix<F, Z>::DenseMatrix(int rows, int cols)
        : boost_dense_matrix(rows, cols),
          usedRanges_(cols, Interval(0, 0)),
          columnBeingEdited_(-1)
    {
        for (int j = 0; j < cols; j++)
        {
            CheckInvariants(j);
        }
    }

    template<typename F, typename Z>
    DenseMatrix<F, Z>::~DenseMatrix()
    {}

    template<typename F, typename Z>
    int
    DenseMatrix<F, Z>::UsedEntries() const
    {
        // use column ranges
        int filledEntries = 0;
        for (int col = 0; col < Columns(); ++col)
        {
            int start, end;
            boost::tie(start, end) = UsedRowRange(col);
            filledEntries += (end - start);
        }
        return filledEntries;
    }

    template<typename F, typename Z>
    int
    DenseMatrix<F, Z>::AllocatedEntries() const
    {
        return Rows() * Columns();
    }

    template<typename F, typename Z>
    void
    DenseMatrix<F, Z>::ToHostMatrix(F** mat, int* rows, int* cols) const
    {
        using boost::numeric::ublas::matrix;
        using boost::numeric::ublas::row_major;

        // TODO(dalexander): make sure SWIG client deallocates this memory -- use %newobject flag
        matrix<Z, row_major> rowMajorPeer(*this);
        *mat = new F[Rows() * Columns()];
        std::copy(rowMajorPeer.data().begin(), rowMajorPeer.data().end(), *mat);
        *rows = Rows();
        *cols = Columns();
    }

    template<typename F, typename Z>
    void
    DenseMatrix<F, Z>::CheckInvariants(int column) const
    {
        // make sure no used entries are outside of the bands
        int start, end;
        boost::tie(start, end) = UsedRowRange(column);
        assert(0 <= start && start <= end && end <= Rows());
        for (int i = 0; i < Rows(); i++)
        {
            if (!(start <= i && i < end))
            {
                assert ((*this)(i, column) == typename boost_dense_matrix::value_type());
            }
        }
    }
}
