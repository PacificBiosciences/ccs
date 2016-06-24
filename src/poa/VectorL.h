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
#include <cassert>
#include <vector>

namespace PacBio {
namespace Consensus {
namespace detail {

// This is necessary because of dumb C++ linker rules.
// see: http://isocpp.org/wiki/faq/templates#template-friends
// and: http://yosefk.com/c++fqa/templates.html#fqa-35.16  (for hilarious
// rejoinder)
template <typename T>
class VectorL;  // fwd
template <typename T>
T Max(const VectorL<T>& v);
template <typename T>
size_t ArgMax(const VectorL<T>& v);

//
// Vector class that stores only a subsequence of the rows
// (SparseVector has a lot of weird functionality built-in, so we can't use it
// without
//  cleaning it up/refactoring it quite a bit)
//
template <typename T>
class VectorL
{
private:
    std::vector<T> storage_;
    size_t beginRow_;
    size_t endRow_;

public:
    VectorL(int beginRow, int endRow, T defaultVal = T())
        : storage_(endRow - beginRow, defaultVal), beginRow_(beginRow), endRow_(endRow)
    {
    }

    T& operator[](size_t pos)
    {
        assert(beginRow_ <= pos && pos < endRow_);
        return storage_[pos - beginRow_];
    }

    const T& operator[](size_t pos) const
    {
        assert(beginRow_ <= pos && pos < endRow_);
        return storage_[pos - beginRow_];
    }

    size_t BeginRow() const { return beginRow_; }
    size_t EndRow() const { return endRow_; }
    friend T Max<>(const VectorL<T>& v);
    friend size_t ArgMax<>(const VectorL<T>& v);
};

using std::max_element;
using std::distance;

template <typename T>
T Max(const VectorL<T>& v)
{
    return *max_element(v.storage_.begin(), v.storage_.end());
}

template <typename T>
size_t ArgMax(const VectorL<T>& v)
{
    return v.beginRow_ +
           distance(v.storage_.begin(), max_element(v.storage_.begin(), v.storage_.end()));
}

}  // namespace detail
}  // namespace Consensus
}  // PacBio
