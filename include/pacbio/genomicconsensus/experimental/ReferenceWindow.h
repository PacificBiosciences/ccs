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

///
/// \brief Overlap
///
/// \return if windows on same reference, and intervals overlap
///
inline bool Overlap(const ReferenceWindow& lhs, const ReferenceWindow& rhs)
{
    if (lhs.name != rhs.name) return false;
    return lhs.interval.Overlaps(rhs.interval);
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
