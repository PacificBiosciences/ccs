// Author: Derek Barnett

#pragma once

#include <algorithm>

#include <pacbio/data/Interval.h>

namespace PacBio {
namespace GenomicConsensus {

inline std::vector<PacBio::Data::Interval> SplitInterval(const PacBio::Data::Interval& source,
                                                         const size_t span, const size_t overhang)
{
    using PacBio::Data::Interval;
    std::vector<Interval> result;

    size_t pos = source.Left();
    while (pos < source.Right()) {
        const auto left = ((pos < overhang) ? 0 : pos - overhang);
        const auto right = pos + span + overhang;
        result.emplace_back(source.Intersect({left, right}));
        pos += span;
    }
    return result;
}

inline std::vector<PacBio::Data::Interval> SplitInterval(const PacBio::Data::Interval& source,
                                                         const size_t span)
{
    return SplitInterval(source, span, 0);
}

}  // namespace GenomicConsensus
}  // namespace PacBio