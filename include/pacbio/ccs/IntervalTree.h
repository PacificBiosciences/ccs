// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
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

// Author: Lance Hepler

#pragma once

#include <algorithm>
#include <set>

#include <pacbio/ccs/Interval.h>

namespace PacBio {
namespace CCS {

class IntervalTree
{
private:
    struct WeakIntervalOrdering
    {
        bool operator()(const Interval& lhs, const Interval& rhs) const
        {
            return lhs.Left() < rhs.Left();
        }
    };

    typedef std::multiset<Interval, WeakIntervalOrdering> StorageType;

public:
    typedef StorageType::iterator       iterator;
    typedef StorageType::const_iterator const_iterator;

    inline
    void Insert(const Interval& interval)
    {
        auto it = storage.insert(interval);

        // check to see if we overlap the previous element,
        // if we do, start our merge loop from there
        if (it != begin())
        {
            const_iterator prev = std::prev(it);

            if (prev->Overlaps(*it))
                it = prev;
        }

        while (it != end())
        {
            const_iterator nx = std::next(it);

            if (nx != end() && it->Overlaps(*nx))
            {
                const Interval u = it->Union(*nx);
                it = storage.erase(it);
                it = storage.erase(it);
                it = storage.insert(it, u);
            }
            else break;
        }
    }

    inline
    IntervalTree Gaps() const
    {
        IntervalTree gaps;

        for (auto it = begin(); it != end(); ++it)
        {
            const_iterator nx = std::next(it);

            if (nx == end())
                break;

            gaps.Insert(Interval(it->Right(), nx->Left()));
        }

        return gaps;
    }

    inline
    IntervalTree Gaps(const Interval& interval) const
    {
        auto left  = begin();
        auto right = end();

        // if we're empty (left == right == end()) or we don't overlap the interval
        // return just the provided interval
        if (left == right ||
            !interval.Overlaps(Interval(left->Left(), (--right)->Right())))
        {
            IntervalTree gaps;
            gaps.Insert(interval);
            return gaps;
        }

        IntervalTree gaps = Gaps();

        if (interval.Left() < left->Left())
            gaps.Insert(Interval(interval.Left(), left->Left()));

        if (right->Right() < interval.Right())
            gaps.Insert(Interval(right->Right(), interval.Right()));

        return gaps;
    }

    inline
    iterator begin()
    {
        return storage.begin();
    }

    inline
    const_iterator begin() const
    {
        return storage.begin();
    }

    inline
    iterator end()
    {
        return storage.end();
    }

    inline
    const_iterator end() const
    {
        return storage.end();
    }

    inline
    size_t size() const
    {
        return storage.size();
    }

private:
    StorageType storage;
};

} // namespace CCS
} // namespace PacBio
