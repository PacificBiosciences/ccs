// Author: Derek Barnett

#pragma once

#include <iostream>
#include <string>
#include <tuple>

#include <pacbio/data/Interval.h>

namespace PacBio {
namespace GenomicConsensus {

struct ReferenceWindow
{
    std::string name;
    PacBio::Data::Interval interval;

    size_t Start(void) const { return interval.Left(); }
    size_t End(void) const { return interval.Right(); }
    size_t Length(void) const { return interval.Length(); }

    bool operator==(const ReferenceWindow& other) const
    {
        return name == other.name && interval == other.interval;
    }

    bool operator!=(const ReferenceWindow& other) const { return !(*this == other); }

    bool operator<(const ReferenceWindow& other) const
    {
        return std::tie(name, interval) < std::tie(other.name, other.interval);
    }
};

inline std::ostream& operator<<(std::ostream& os, const ReferenceWindow& window)
{
    os << window.name << ' ' << window.interval;
    return os;
}

}  // namespace GenomicConsensus
}  // namespace PacBio
