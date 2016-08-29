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

#pragma once

#include <cstddef>
#include "pacbio/consensus/AbstractMatrix.h"

namespace PacBio {
namespace Consensus {

//
// BasicDenseMatrix is a *basic* dense matrix, for use as an
// intermediate in matrix viewing operations (not in production code).
//
// It does not fully implement the interface that would be required to
// drop it in as a replacement for ScaledMatrix in the production code
// (i.e. for the recursor, etc.).  ConsensusCore did offer such a
// matrix, "DenseMatrix", and we could consider resurrecting such a
// class in the future.
//
class BasicDenseMatrix : public AbstractMatrix
{
private:
    size_t nCols_;
    size_t nRows_;
    double* entries_;

public:  // structors
    BasicDenseMatrix(size_t rows, size_t cols);
    BasicDenseMatrix(const BasicDenseMatrix& other) = delete;
    ~BasicDenseMatrix();

public:  // Size information
    size_t Rows() const;
    size_t Columns() const;

public:  // accessors
    double& operator()(size_t i, size_t j) const;

public:  // AbstractMatrix interface
    void ToHostMatrix(double** mat, int* rows, int* cols) const;
    size_t UsedEntries() const;
    float UsedEntriesRatio() const;
    size_t AllocatedEntries() const;
};

inline size_t BasicDenseMatrix::Rows() const { return nRows_; }

inline size_t BasicDenseMatrix::Columns() const { return nCols_; }

inline double& BasicDenseMatrix::operator()(size_t i, size_t j) const
{
    return entries_[i * Columns() + j];
}
}
}
