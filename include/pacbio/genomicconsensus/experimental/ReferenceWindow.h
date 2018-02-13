// Copyright (c) 2017, Pacific Biosciences of California, Inc.
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

// Author: Derek Barnett

#pragma once

#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include <pacbio/data/Interval.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

///
/// \brief The ReferenceWindow struct
///
struct ReferenceWindow
{
    std::string name;
    PacBio::Data::Interval interval;

    // helpers for interval access
    size_t Start(void) const { return interval.Left(); }
    size_t End(void) const { return interval.Right(); }
    size_t Length(void) const { return interval.Length(); }
};

///
/// \brief AreContiguous
/// \param windows
/// \return
///
inline bool AreContiguous(const std::vector<ReferenceWindow>& windows)
{
    //
    // Predicate that determines whether the reference/scaffold windows
    // are contiguous.
    //

    std::string lastName;
    size_t lastEnd = 0;
    for (const auto& win : windows) {
        if ((!lastName.empty() && win.name != lastName) ||
            (lastEnd != 0 && win.Start() != lastEnd)) {
            return false;
        }
        lastEnd = win.End();
        lastName = win.name;
    }
    return true;
}

///
/// \brief AreContiguous
///
/// Helper for checking 2 windows.
///
/// \param window1
/// \param window2
/// \return
///
inline bool AreContiguous(const ReferenceWindow& lhs, const ReferenceWindow& rhs)
{
    return AreContiguous({lhs, rhs});
}

inline bool operator==(const ReferenceWindow& lhs, const ReferenceWindow& rhs)
{
    return std::tie(lhs.name, lhs.interval) == std::tie(rhs.name, rhs.interval);
}

inline bool operator!=(const ReferenceWindow& lhs, const ReferenceWindow& rhs)
{
    return !(lhs == rhs);
}

inline bool operator<(const ReferenceWindow& lhs, const ReferenceWindow& rhs)
{
    return std::tie(lhs.name, lhs.interval) < std::tie(rhs.name, rhs.interval);
}

inline std::ostream& operator<<(std::ostream& os, const ReferenceWindow& window)
{
    os << window.name << ' ' << window.interval;
    return os;
}

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
