// Copyright (c) 2011-2016, Pacific Biosciences of California, Inc.
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

#include "SparseMatrix.h"

namespace PacBio {
namespace Consensus {

SparseMatrix::SparseMatrix(const size_t rows, const size_t cols)
    : columns_(cols)
    , nCols_(cols)
    , nRows_(rows)
    , columnBeingEdited_(std::numeric_limits<size_t>::max())
    , usedRanges_(cols, std::make_pair(0, 0))
{
}

SparseMatrix::SparseMatrix(const SparseMatrix& other)
    : columns_(other.nCols_)
    , nCols_(other.nCols_)
    , nRows_(other.nRows_)
    , columnBeingEdited_(other.columnBeingEdited_)
    , usedRanges_(other.usedRanges_)
{
    for (size_t j = 0; j < nCols_; j++)
        if (other.columns_[j])
            columns_[j] = std::unique_ptr<SparseVector>(new SparseVector(*other.columns_[j]));
}

SparseMatrix::~SparseMatrix() {}
void SparseMatrix::Reset(const size_t rows, const size_t cols)
{
    std::vector<std::unique_ptr<SparseVector>>(cols).swap(columns_);
    nCols_ = cols;
    nRows_ = rows;
    std::vector<std::pair<size_t, size_t>>(cols, std::make_pair(0, 0)).swap(usedRanges_);
    columnBeingEdited_ = std::numeric_limits<size_t>::max();
}

size_t SparseMatrix::UsedEntries() const
{
    // use column ranges
    size_t filledEntries = 0;
    for (size_t col = 0; col < Columns(); ++col) {
        size_t start, end;
        std::tie(start, end) = UsedRowRange(col);
        filledEntries += (end - start);
    }
    return filledEntries;
}

size_t SparseMatrix::AllocatedEntries() const
{
    size_t sum = 0;
    for (size_t j = 0; j < nCols_; j++) {
        sum += (columns_[j] != NULL ? columns_[j]->AllocatedEntries() : 0);
    }
    return sum;
}

void SparseMatrix::ToHostMatrix(double** mat, size_t* rows, size_t* cols) const
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    *mat = new double[Rows() * Columns()];
    *rows = Rows();
    *cols = Columns();
    for (size_t i = 0; i < Rows(); i++) {
        for (size_t j = 0; j < Columns(); j++) {
            (*mat)[i * Columns() + j] = IsAllocated(i, j) ? Get(i, j) : nan;
        }
    }
}

void SparseMatrix::CheckInvariants(size_t column) const
{
#ifndef NDEBUG
    for (size_t j = 0; j < nCols_; j++) {
        if (columns_[j] != NULL) columns_[j]->CheckInvariants();
    }
#endif
}

}  // namespace Consensus
}  // namespace PacBio
