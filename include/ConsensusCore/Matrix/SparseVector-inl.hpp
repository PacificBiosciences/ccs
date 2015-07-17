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

#include <algorithm>
#include <cassert>
#include <vector>

#define PADDING          8
#define SHRINK_THRESHOLD 0.8

namespace ConsensusCore
{
    template<typename F, typename Z>
    inline
    SparseVector<F, Z>::SparseVector(int logicalLength, int beginRow, int endRow)
        : logicalLength_(logicalLength)
        , allocatedBeginRow_(std::max(beginRow - PADDING, 0))
        , allocatedEndRow_(std::min(endRow + PADDING, logicalLength_))
        , storage_(allocatedEndRow_ - allocatedBeginRow_, Z())
        , nReallocs_(0)
    {
        assert(beginRow >= 0      &&
               beginRow <= endRow &&
               endRow   <= logicalLength);
        DEBUG_ONLY(CheckInvariants());
    }

    template<typename F, typename Z>
    inline
    SparseVector<F, Z>::SparseVector(const SparseVector<F, Z>& other)
        : logicalLength_(other.logicalLength_)
        , allocatedBeginRow_(other.allocatedBeginRow_)
        , allocatedEndRow_(other.allocatedEndRow_)
        , storage_(other.storage_)
        , nReallocs_(0)
    {
        DEBUG_ONLY(CheckInvariants());
    }

    template<typename F, typename Z>
    inline void
    SparseVector<F, Z>::ResetForRange(int beginRow, int endRow)
    {
        // Allows reuse.  Destructive.
        DEBUG_ONLY(CheckInvariants());
        assert(beginRow >= 0      &&
               beginRow <= endRow &&
               endRow   <= logicalLength_);
        int newAllocatedBegin =  std::max(beginRow - PADDING, 0);
        int newAllocatedEnd   =  std::min(endRow   + PADDING, logicalLength_);
        if ((newAllocatedEnd - newAllocatedBegin) > (allocatedEndRow_ - allocatedBeginRow_))
        {
            storage_.resize(newAllocatedEnd - newAllocatedBegin);
            nReallocs_++;
            Clear();
        }
        else if ((newAllocatedEnd - newAllocatedBegin) <
                 static_cast<int>(SHRINK_THRESHOLD * (allocatedEndRow_ - allocatedBeginRow_)))
        {
            // use swap trick to free allocated but unused memory,
            // see: http://stackoverflow.com/questions/253157/how-to-downsize-stdvector
            std::vector<F>(newAllocatedEnd - newAllocatedBegin, Z()).swap(storage_);
            nReallocs_++;
        }
        else
        {
            Clear();
        }
        allocatedBeginRow_ = newAllocatedBegin;
        allocatedEndRow_   = newAllocatedEnd;
        DEBUG_ONLY(CheckInvariants());
    }

    template<typename F, typename Z>
    inline void
    SparseVector<F, Z>::ExpandAllocated(int newAllocatedBegin, int newAllocatedEnd)
    {
        // Expands allocated storage while preserving the contents.
        DEBUG_ONLY(CheckInvariants());
        assert(newAllocatedBegin >= 0                  &&
               newAllocatedBegin <= newAllocatedEnd    &&
               newAllocatedEnd   <= logicalLength_);
        assert(newAllocatedBegin <= allocatedBeginRow_ &&
               newAllocatedEnd   >= allocatedEndRow_);
        // Resize the underlying storage.
        storage_.resize(newAllocatedEnd - newAllocatedBegin);
        // Use memmove to robustly relocate the old data (handles overlapping ranges).
        //   Data is at:
        //      storage[0 ... (end - begin) )
        //   Must be moved to:
        //      storage[(begin - newBegin) ... (end - newBegin)]
        memmove(&storage_[allocatedBeginRow_ - newAllocatedBegin],
                &storage_[0],
                (allocatedEndRow_ - allocatedBeginRow_) * sizeof(F)); // NOLINT
        // "Zero"-fill the allocated but unused space.
        std::fill(storage_.begin(),
                  storage_.begin() + (allocatedBeginRow_ - newAllocatedBegin),
                  Z());
        std::fill(storage_.begin() + (allocatedEndRow_- newAllocatedBegin),
                  storage_.end(),
                  Z());
        // Update pointers.
        allocatedBeginRow_ = newAllocatedBegin;
        allocatedEndRow_   = newAllocatedEnd;
        nReallocs_++;
        DEBUG_ONLY(CheckInvariants());
    }

    template<typename F, typename Z>
    inline bool
    SparseVector<F, Z>::IsAllocated(int i) const
    {
        assert(i >= 0 && i < logicalLength_);
        return i >= allocatedBeginRow_ && i < allocatedEndRow_;
    }

    template<typename F, typename Z>
    inline const F&
    SparseVector<F, Z>::operator()(int i) const
    {
        if (IsAllocated(i))
        {
            return storage_[i - allocatedBeginRow_];
        }
        else
        {
            static const F emptyCell_ = Z();
            return emptyCell_;
        }
    }

    template<typename F, typename Z>
    inline F
    SparseVector<F, Z>::Get(int i) const
    {
        return (*this)(i);
    }

    template<typename F, typename Z>
    inline void
    SparseVector<F, Z>::Set(int i, F v)
    {
        using std::max;
        using std::min;

        DEBUG_ONLY(CheckInvariants());
        assert (i >= 0 && i < logicalLength_);
        if (!IsAllocated(i))
        {
            int newBeginRow = max(min(i - PADDING, allocatedBeginRow_), 0);
            int newEndRow   = min(max(i + PADDING, allocatedEndRow_), logicalLength_);
            ExpandAllocated(newBeginRow, newEndRow);
        }
        storage_[i - allocatedBeginRow_] = v;
        DEBUG_ONLY(CheckInvariants());
    }

    template<>
    inline __m128
    SparseVector<float, lvalue<float>>::Get4(int i) const
    {
        assert(i >= 0 && i < logicalLength_ - 3);
        if (i >= allocatedBeginRow_ && i < allocatedEndRow_ - 3)
        {
            return _mm_loadu_ps(&storage_[i-allocatedBeginRow_]);
        }
        else
        {
            return _mm_set_ps(Get(i+3), Get(i+2), Get(i+1), Get(i+0));
        }
    }

    template<typename F, typename Z>
    inline __m128
    SparseVector<F, Z>::Get4(int i) const
    {
        throw std::runtime_error("cannot perform Get4 with non-f32 type");
    }

    template<>
    inline void
    SparseVector<float, lvalue<float>>::Set4(int i, __m128 v4)
    {
        assert(i >= 0 && i < logicalLength_ - 3);
        if (i >= allocatedBeginRow_ && i < allocatedEndRow_ - 3)
        {
            _mm_storeu_ps(&storage_[i-allocatedBeginRow_], v4);
        }
        else
        {
            float vbuf[4];
            _mm_storeu_ps(vbuf, v4);
            Set(i+0, vbuf[0]);
            Set(i+1, vbuf[1]);
            Set(i+2, vbuf[2]);
            Set(i+3, vbuf[3]);
        }
    }

    template<typename F, typename Z>
    inline void
    SparseVector<F, Z>::Set4(int i, __m128 v4)
    {
        throw std::runtime_error("cannot perform Set4 with non-f32 type");
    }

    template<typename F, typename Z>
    inline void
    SparseVector<F, Z>::Clear()
    {
        std::fill(storage_.begin(), storage_.end(), Z());
    }

    template<typename F, typename Z>
    inline int
    SparseVector<F, Z>::AllocatedEntries() const
    {
        // We want the real memory usage.  std::vector is holding some memory back
        // from us.
        return storage_.capacity();
    }

    template<typename F, typename Z>
    inline void
    SparseVector<F, Z>::CheckInvariants() const
    {
        assert(logicalLength_ >= 0);
        assert(0 <= allocatedBeginRow_ && allocatedBeginRow_ < logicalLength_);
        assert(0 <= allocatedEndRow_ && allocatedEndRow_ <= logicalLength_);
        assert(allocatedBeginRow_ <= allocatedEndRow_);
        assert((allocatedEndRow_ - allocatedBeginRow_) <= (signed)storage_.size());
    }
}

#undef PADDING
#undef SHRINK_THRESHOLD
