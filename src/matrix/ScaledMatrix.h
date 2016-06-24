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

#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

#include "SparseMatrix.h"

namespace PacBio {
namespace Consensus {

class ScaledMatrix : public SparseMatrix
{
public:
    enum Direction
    {
        FORWARD,
        REVERSE
    };

public:  // constructor/destructor
    ScaledMatrix(int rows, int cols, Direction dir);
    ScaledMatrix(const ScaledMatrix& other);
    ~ScaledMatrix(void) = default;

public:
    void Reset(size_t rows, size_t cols);
    Direction SetDirection(Direction dir);

public:  // nullability
    static const ScaledMatrix& Null();

public:  // information about entries filled by column
    template <bool maxProvided>
    void FinishEditingColumn(int j, int usedBegin, int usedEnd, double max_val = 0.0);

public:  // Scaling and normalization
    double GetLogScale(int j) const;
    double GetLogProdScales(int s, int e) const;
    double GetLogProdScales() const;

private:
    std::vector<double> logScalars_;
    Direction dir_;
};

std::ostream& operator<<(std::ostream&, const ScaledMatrix&);

inline const ScaledMatrix& ScaledMatrix::Null()
{
    static ScaledMatrix* nullObj = new ScaledMatrix(0, 0, FORWARD);
    return *nullObj;
}

template <bool maxProvided>
inline void ScaledMatrix::FinishEditingColumn(const int j, const int usedBegin, const int usedEnd,
                                              double max_val)
{
    // get the constant to scale by
    if (!maxProvided) {
        max_val = 0.0;
        for (int i = usedBegin; i < usedEnd; ++i) {
            max_val = std::max(max_val, SparseMatrix::Get(i, j));
        }
    }

    // cumsum stuff
    double last = 0.0;
    if (dir_ == FORWARD && j > 0)
        last = logScalars_[j - 1];
    else if (dir_ == REVERSE && j + 1 < Columns())
        last = logScalars_[j + 1];

    // set it
    if (max_val != 0.0 && max_val != 1.0) {
        for (int i = usedBegin; i < usedEnd; ++i) {
            SparseMatrix::Set(i, j, SparseMatrix::Get(i, j) / max_val);
        }
        logScalars_[j] = last + std::log(max_val);
    } else {
        logScalars_[j] = last;
    }

    SparseMatrix::FinishEditingColumn(j, usedBegin, usedEnd);
}

inline double ScaledMatrix::GetLogScale(int j) const { return logScalars_[j]; }
inline double ScaledMatrix::GetLogProdScales(int beginColumn, int endColumn) const
{
    double f, l;
    if (dir_ == FORWARD) {
        f = (beginColumn > 0) ? logScalars_[beginColumn - 1] : 0.0;
        l = (endColumn > 0) ? logScalars_[endColumn - 1] : 0.0;
        return l - f;
    }  // else dir_ == REVERSE
    f = logScalars_[beginColumn];
    l = (endColumn < Columns()) ? logScalars_[endColumn] : 0.0;
    return f - l;
}

inline double ScaledMatrix::GetLogProdScales() const
{
    if (dir_ == FORWARD) return logScalars_.back();
    return logScalars_.front();
}

}  // namespace Consensus
}  // namespace PacBio
