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

#include <algorithm>
#include <cassert>
#include <cstring>
#include <utility>
#include <vector>

namespace PacBio {
namespace Consensus {

class SparseVector
{
public:  // Constructor, destructor
    SparseVector(size_t logicalLength, size_t beginRow, size_t endRow);
    SparseVector(const SparseVector& other);
    ~SparseVector() = default;

    // Ensures there is enough allocated storage to
    // hold entries for at least [beginRow, endRow) (plus padding);
    // clears existing entries.
    void ResetForRange(size_t beginRow, size_t endRow);

public:
    const double& operator()(size_t i) const;
    bool IsAllocated(size_t i) const;
    double Get(size_t i) const;
    void Set(size_t i, double v);
    void Clear();

public:
    size_t AllocatedEntries() const;
    void CheckInvariants() const;

private:
    // Expand the range of rows for which we have backing storage,
    // while preserving contents.  The arguments will become the
    // new allocated bounds, so caller should add padding if desired
    // before calling.
    void ExpandAllocated(size_t newAllocatedBegin, size_t newAllocatedEnd);

private:
    // the "logical" length of the vector, of which only
    // a subset of entries are actually allocated
    size_t logicalLength_;

    // row numbers in the abstraction we are presenting
    size_t allocatedBeginRow_;
    size_t allocatedEndRow_;

    // the storage
    std::vector<double> storage_;

    // analytics
    size_t nReallocs_;
};

namespace {  // anonymous

constexpr size_t PADDING = 8;
constexpr double SHRINK_THRESHOLD = 0.8;

}  // anonymous namespace

inline SparseVector::SparseVector(size_t logicalLength, size_t beginRow, size_t endRow)
    : logicalLength_(logicalLength)
    , allocatedBeginRow_((beginRow > PADDING) ? beginRow - PADDING : 0)
    , allocatedEndRow_(std::min(endRow + PADDING, logicalLength_))
    , storage_(allocatedEndRow_ - allocatedBeginRow_, 0.0)
    , nReallocs_(0)
{
    assert(beginRow >= 0 && beginRow <= endRow && endRow <= logicalLength);
    CheckInvariants();
}

inline SparseVector::SparseVector(const SparseVector& other)
    : logicalLength_(other.logicalLength_)
    , allocatedBeginRow_(other.allocatedBeginRow_)
    , allocatedEndRow_(other.allocatedEndRow_)
    , storage_(other.storage_)
    , nReallocs_(0)
{
    CheckInvariants();
}

inline void SparseVector::ResetForRange(size_t beginRow, size_t endRow)
{
    // Allows reuse.  Destructive.
    CheckInvariants();
    assert(beginRow >= 0 && beginRow <= endRow && endRow <= logicalLength_);
    size_t newAllocatedBegin = (beginRow > PADDING) ? beginRow - PADDING : 0;
    size_t newAllocatedEnd = std::min(endRow + PADDING, logicalLength_);
    if ((newAllocatedEnd - newAllocatedBegin) > (allocatedEndRow_ - allocatedBeginRow_)) {
        storage_.resize(newAllocatedEnd - newAllocatedBegin);
        nReallocs_++;
        Clear();
    } else if ((newAllocatedEnd - newAllocatedBegin) <
               static_cast<size_t>(SHRINK_THRESHOLD * (allocatedEndRow_ - allocatedBeginRow_))) {
        // use swap trick to free allocated but unused memory,
        // see:
        // http://stackoverflow.com/questions/253157/how-to-downsize-stdvector
        std::vector<double>(newAllocatedEnd - newAllocatedBegin, 0.0).swap(storage_);
        nReallocs_++;
    } else {
        Clear();
    }
    allocatedBeginRow_ = newAllocatedBegin;
    allocatedEndRow_ = newAllocatedEnd;
    CheckInvariants();
}

inline void SparseVector::ExpandAllocated(size_t newAllocatedBegin, size_t newAllocatedEnd)
{
    // Expands allocated storage while preserving the contents.
    CheckInvariants();
    assert(newAllocatedBegin >= 0 && newAllocatedBegin <= newAllocatedEnd &&
           newAllocatedEnd <= logicalLength_);
    assert(newAllocatedBegin <= allocatedBeginRow_ && newAllocatedEnd >= allocatedEndRow_);
    // Resize the underlying storage.
    storage_.resize(newAllocatedEnd - newAllocatedBegin);
    // Use memmove to robustly relocate the old data (handles overlapping
    // ranges).
    //   Data is at:
    //      storage[0 ... (end - begin) )
    //   Must be moved to:
    //      storage[(begin - newBegin) ... (end - newBegin)]
    memmove(&storage_[allocatedBeginRow_ - newAllocatedBegin], &storage_[0],
            (allocatedEndRow_ - allocatedBeginRow_) * sizeof(double));  // NOLINT
    // "Zero"-fill the allocated but unused space.
    std::fill(storage_.begin(), storage_.begin() + (allocatedBeginRow_ - newAllocatedBegin), 0.0);
    std::fill(storage_.begin() + (allocatedEndRow_ - newAllocatedBegin), storage_.end(), 0.0);
    // Update pointers.
    allocatedBeginRow_ = newAllocatedBegin;
    allocatedEndRow_ = newAllocatedEnd;
    nReallocs_++;
    CheckInvariants();
}

inline bool SparseVector::IsAllocated(size_t i) const
{
    assert(i >= 0 && i < logicalLength_);
    return i >= allocatedBeginRow_ && i < allocatedEndRow_;
}

inline const double& SparseVector::operator()(size_t i) const
{
    if (IsAllocated(i)) {
        return storage_[i - allocatedBeginRow_];
    } else {
        static const double emptyCell_ = 0.0;
        return emptyCell_;
    }
}

inline double SparseVector::Get(size_t i) const { return (*this)(i); }
inline void SparseVector::Set(size_t i, double v)
{
    using std::max;
    using std::min;

    CheckInvariants();
    assert(i >= 0 && i < logicalLength_);
    if (!IsAllocated(i)) {
        size_t newBeginRow = min((i > PADDING) ? i - PADDING : 0, allocatedBeginRow_);
        size_t newEndRow = min(max(i + PADDING, allocatedEndRow_), logicalLength_);
        ExpandAllocated(newBeginRow, newEndRow);
    }
    storage_[i - allocatedBeginRow_] = v;
    CheckInvariants();
}

inline void SparseVector::Clear() { std::fill(storage_.begin(), storage_.end(), 0.0); }
inline size_t SparseVector::AllocatedEntries() const
{
    // We want the real memory usage.  std::vector is holding some memory back
    // from us.
    return storage_.capacity();
}

inline void SparseVector::CheckInvariants() const
{
    assert(logicalLength_ >= 0);
    assert(0 <= allocatedBeginRow_ && allocatedBeginRow_ < logicalLength_);
    assert(0 <= allocatedEndRow_ && allocatedEndRow_ <= logicalLength_);
    assert(allocatedBeginRow_ <= allocatedEndRow_);
    assert((allocatedEndRow_ - allocatedBeginRow_) <= storage_.size());
}

}  // namespace Consensus
}  // namespace PacBio
