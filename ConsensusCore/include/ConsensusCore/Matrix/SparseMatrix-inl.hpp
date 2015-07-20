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

#include <algorithm>
#include <boost/tuple/tuple.hpp>
#include <cassert>

#include <ConsensusCore/Interval.hpp>
#include <ConsensusCore/LValue.hpp>

namespace ConsensusCore {
    //
    // Nullability
    //
    template<typename F, typename Z>
    inline const SparseMatrix<F, Z>&
    SparseMatrix<F, Z>::Null()
    {
        static SparseMatrix<F, Z>* nullObj = new SparseMatrix<F, Z>(0, 0);
        return *nullObj;
    }

    template<typename F, typename Z>
    inline bool
    SparseMatrix<F, Z>::IsNull() const
    {
        return (Rows() == 0 && Columns() == 0);
    }

    //
    // Size information
    //
    template<typename F, typename Z>
    inline const int
    SparseMatrix<F, Z>::Rows() const
    {
        return nRows_;
    }

    template<typename F, typename Z>
    inline const int
    SparseMatrix<F, Z>::Columns() const
    {
        return nCols_;
    }

    //
    // Entry range queries per column
    //
    template<typename F, typename Z>
    inline void
    SparseMatrix<F, Z>::StartEditingColumn(int j, int hintBegin, int hintEnd)
    {
        assert(columnBeingEdited_ == -1);
        columnBeingEdited_ = j;
        if (columns_[j] != NULL)
        {
            columns_[j]->ResetForRange(hintBegin, hintEnd);
        } else {
            columns_[j] = new SparseVector<F, Z>(Rows(), hintBegin, hintEnd);
        }
    }

    template<typename F, typename Z>
    inline void
    SparseMatrix<F, Z>::FinishEditingColumn(int j, int usedRowsBegin, int usedRowsEnd)
    {
        assert(columnBeingEdited_ == j);
        usedRanges_[j] = Interval(usedRowsBegin, usedRowsEnd);
        DEBUG_ONLY(CheckInvariants(columnBeingEdited_));
        columnBeingEdited_ = -1;
    }

    template<typename F, typename Z>
    inline Interval
    SparseMatrix<F, Z>::UsedRowRange(int j) const
    {
        assert(0 <= j && j < (int)usedRanges_.size());
        return usedRanges_[j];
    }

    template<typename F, typename Z>
    inline bool
    SparseMatrix<F, Z>::IsColumnEmpty(int j) const
    {
        assert(0 <= j && j < (int)usedRanges_.size());
        return (usedRanges_[j].Begin >= usedRanges_[j].End);
    }

    //
    // Accessors
    //
    template<typename F, typename Z>
    inline const F&
    SparseMatrix<F, Z>::operator() (int i, int j) const
    {
        if (columns_[j] == NULL)
        {
            static const F emptyCell = Z();
            return emptyCell;
        }
        else
        {
            return (*columns_[j])(i);
        }
    }

    template<typename F, typename Z>
    inline bool
    SparseMatrix<F, Z>::IsAllocated(int i, int j) const
    {
        return columns_[j] != NULL && columns_[j]->IsAllocated(i);
    }

    template<typename F, typename Z>
    inline F
    SparseMatrix<F, Z>::Get(int i, int j) const
    {
        return (*this)(i, j);
    }

    template<typename F, typename Z>
    inline void
    SparseMatrix<F, Z>::Set(int i, int j, F v)
    {
        assert(columnBeingEdited_ == j);
        columns_[j]->Set(i, v);
    }

    template<typename F, typename Z>
    inline void
    SparseMatrix<F, Z>::ClearColumn(int j)
    {
        usedRanges_[j] = Interval(0, 0);
        columns_[j]->Clear();
        DEBUG_ONLY(CheckInvariants(j);)
    }

    //
    // SSE
    //
    template<>
    inline __m128
    SparseMatrix<float, lvalue<float>>::Get4(int i, int j) const
    {
        if (columns_[j] == NULL)
        {
            return Zero4<lvalue<float>>();
        }
        else
        {
            return columns_[j]->Get4(i);
        }
    }

    template<typename F, typename Z>
    inline __m128
    SparseMatrix<F, Z>::Get4(int i, int j) const
    {
        throw std::runtime_error("cannot perform Get4 with non-f32 type");
    }

    template<>
    inline void
    SparseMatrix<float, lvalue<float>>::Set4(int i, int j, __m128 v4)
    {
        assert(columnBeingEdited_ == j);
        columns_[j]->Set4(i, v4);
    }

    template<typename F, typename Z>
    inline void
    SparseMatrix<F, Z>::Set4(int i, int j, __m128 v4)
    {
        throw std::runtime_error("cannot perform Get4 with non-f32 type");
    }

    //
    // Performance insensitive routines are not inlined
    //
    template<typename F, typename Z>
    SparseMatrix<F, Z>::SparseMatrix(int rows, int cols)
        : columns_(cols),  nCols_(cols), nRows_(rows), columnBeingEdited_(-1),
          usedRanges_(cols, Interval(0, 0))
    {
        for (int j = 0; j < nCols_; j++)
        {
            columns_[j] = NULL;
        }
    }

    template<typename F, typename Z>
    SparseMatrix<F, Z>::SparseMatrix(const SparseMatrix& other)
        : columns_(other.nCols_),
          nCols_(other.nCols_),
          nRows_(other.nRows_),
          columnBeingEdited_(other.columnBeingEdited_),
          usedRanges_(other.usedRanges_)
    {
        for (int j = 0; j < nCols_; j++)
        {
            if (other.columns_[j] != NULL)
            {
                columns_[j] = new SparseVector<F, Z>(*other.columns_[j]);
            }
            else
            {
                columns_[j] = NULL;
            }
        }
    }

    template<typename F, typename Z>
    SparseMatrix<F, Z>::~SparseMatrix()
    {
        for (int j = 0; j < nCols_; j++)
        {
            if (columns_[j] != NULL) delete columns_[j];
        }
    }

    template<typename F, typename Z>
    int
    SparseMatrix<F, Z>::UsedEntries() const
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
    SparseMatrix<F, Z>::AllocatedEntries() const
    {
        int sum = 0;
        for (int j = 0; j < nCols_; j++)
        {
            sum += (columns_[j] != NULL ?
                    columns_[j]->AllocatedEntries() : 0);
        }
        return sum;
    }

    template<typename F, typename Z>
    void
    SparseMatrix<F, Z>::ToHostMatrix(F** mat, int* rows, int* cols) const
    {
        const F nan = std::numeric_limits<F>::quiet_NaN();
        *mat = new F[Rows() * Columns()];
        *rows = Rows();
        *cols = Columns();
        for (int i = 0; i < Rows(); i++) {
            for (int j = 0; j < Columns(); j++) {
                (*mat)[i * Columns() + j] = IsAllocated(i, j) ? Get(i, j) : nan;
            }
        }
    }

    template<typename F, typename Z>
    void
    SparseMatrix<F, Z>::CheckInvariants(int column) const
    {
        for (int j = 0; j < nCols_; j++)
         {
             if (columns_[j] != NULL) columns_[j]->CheckInvariants();
         }
    }
}
