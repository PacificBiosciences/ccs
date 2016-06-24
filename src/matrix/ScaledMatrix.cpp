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

#include <algorithm>
#include <cmath>
#include <iomanip>

#include "ScaledMatrix.h"

namespace PacBio {
namespace Consensus {

ScaledMatrix::ScaledMatrix(int rows, int cols, Direction dir)
    : SparseMatrix(rows, cols), logScalars_(cols, 0.0), dir_{dir}
{
}
ScaledMatrix::ScaledMatrix(const ScaledMatrix& other)
    : SparseMatrix(other), logScalars_(other.logScalars_)
{
}

void ScaledMatrix::Reset(size_t rows, size_t cols)
{
    std::vector<double>(cols, 0.0).swap(logScalars_);
    SparseMatrix::Reset(rows, cols);
}

ScaledMatrix::Direction ScaledMatrix::SetDirection(const Direction dir)
{
    const Direction res = dir_;
    dir_ = dir;
    std::fill(logScalars_.begin(), logScalars_.end(), 0.0);
    return res;
}

std::ostream& operator<<(std::ostream& os, const ScaledMatrix& mat)
{
    os << "MATRIX (" << mat.Rows() << ", " << mat.Columns() << ") BEGIN" << std::endl;
    os << std::fixed << std::setprecision(4);
    for (size_t i = 0; i < mat.Rows(); ++i) {
        os << ' ';
        for (size_t j = 0; j < mat.Columns(); ++j) {
            os << ' ' << std::setw(9) << std::log(mat.Get(i, j)) + mat.GetLogScale(j);
        }
        os << std::endl;
    }
    os << "END" << std::endl;
    return os;
}

}  // namespace Consensus
}  // namespace PacBio
