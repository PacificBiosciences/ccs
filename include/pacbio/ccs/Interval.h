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
#include <cassert>
#include <functional>
#include <iosfwd>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

namespace PacBio {
namespace CCS {

class Interval
{
public:
    inline
    Interval()
        : left{0}
        , right{0}
    {
        Validate();
    }

    inline
    Interval(const size_t l, const size_t r)
        : left{l}
        , right{r}
    {
        Validate();
    }

    inline
    void Reset(const size_t l, const size_t r)
    {
        left = l;
        right = r;
        Validate();
    }

    inline
    void Reset(const Interval& other)
    {
        left = other.Left();
        right = other.Right();
        Validate();
    }

    inline
    size_t Length() const
    {
        return right - left;
    }

    inline
    size_t Left() const
    {
        return left;
    }

    inline
    size_t Right() const
    {
        return right;
    }

    // this also handles adjacency
    inline
    bool Overlaps(const Interval& other) const
    {
        // if the left of one is in the range of the other
        return (other.Left() <= left && left <= other.Right()) ||
               (left <= other.Left() && other.Left() <= right);
    }

    inline
    bool Contains(const size_t value) const
    {
        return (left <= value) && (value < right);
    }

    inline
    bool Covers(const Interval& other) const
    {
        return (Intersect(other) == other);
    }

    inline
    Interval Intersect(const Interval& other) const
    {
        if (!Overlaps(other))
            throw std::invalid_argument("interval to intersect does not overlap!");

        return Interval(std::max(left, other.Left()), std::min(right, other.Right()));
    }

    inline
    Interval Union(const Interval& other) const
    {
        if (!Overlaps(other))
            throw std::invalid_argument("interval to merge does not overlap!");

        return Interval(std::min(left, other.Left()), std::max(right, other.Right()));
    }

    inline
    bool operator<(const Interval& other) const
    {
        if (left < other.Left())
            return true;

        if (left == other.Left())
            return right < other.Right();

        return false;
    }

    inline
    bool operator==(const Interval& other) const
    {
        return (left == other.Left() && right == other.Right());
    }

    inline
    bool operator!=(const Interval& other) const
    {
        return !(*this == other);
    }

    inline
    Interval& operator++()
    {
        ++left;
        return *this;
    }

    inline
    Interval operator++(int)
    {
        return Interval(left++, right);
    }

    inline
    size_t operator*() const
    {
        return left;
    }

    inline
    Interval begin() const
    {
        return Interval(left, right);
    }

    inline
    Interval end() const
    {
        return Interval(right, right);
    }

    inline
    operator std::tuple<size_t, size_t>() const
    {
        return std::make_tuple(left, right);
    }

    inline
    operator std::tuple<size_t&, size_t&>()
    {
        return std::make_tuple(std::ref(left), std::ref(right));
    }

    operator std::string() const;
    friend std::ostream& operator<<(std::ostream& os, const Interval& interval);

    static
    Interval FromString(const std::string& str)
    {
        try
        {
            std::vector<std::string> components;
            boost::split(components, str, boost::is_any_of("-"));
            if (components.size() == 1)
            {
                size_t left = boost::lexical_cast<size_t>(components[0]);
                return Interval(left, left + 1);
            }
            else if (components.size() == 2)
            {
                size_t left  = boost::lexical_cast<size_t>(components[0]);
                size_t right = boost::lexical_cast<size_t>(components[1]);
                // if right <= left, we have an invalid interval, fall through
                if (left < right)
                    return Interval(left, right + 1);
            }
        }
        catch (...)
        { }
        throw std::invalid_argument("invalid Interval specification");
    }

private:
    size_t left;
    size_t right;

    inline
    void Validate() const
    {
        assert(left <= right);
    }
};

} // namespace CCS
} // namespace PacBio
