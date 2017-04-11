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

#include <stdio.h>

namespace PacBio {
namespace Consensus {

/// AbstractMatrix is a superclass of the matrix types used in the arrow
/// banded dynamic programming.  It exposes a minimal interface only
/// intended for diagnostic purposes (looking at a matrix from Python,
/// seeing how well the banding is working, ...).  No matrix implementation
/// details are exposed---one can think of this as effectively an opaque
/// data type.
class AbstractMatrix
{
public:
    /// Method SWIG clients can use to get a native matrix (e.g. Numpy)
    /// mat must be filled as a ROW major matrix.
    /// N.B.: Needs int, not size_t dimensions, for SWIG/numpy
    virtual void ToHostMatrix(double** mat, int* rows, int* cols) const = 0;

public:
    // Methods for inquiring about matrix occupancy.
    virtual size_t UsedEntries() const = 0;
    virtual float UsedEntriesRatio() const = 0;
    virtual size_t AllocatedEntries() const = 0;  // an entry may be allocated but not used

public:
    virtual ~AbstractMatrix() {}
};

}  // namespace Consensus
}  // namespace PacBio
