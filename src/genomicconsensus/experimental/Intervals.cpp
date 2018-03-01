// Author: Derek Barnett

#include <pacbio/genomicconsensus/experimental/Intervals.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>

#include <pacbio/consensus/Coverage.h>
#include <pbcopper/logging/Logging.h>

#include <pacbio/genomicconsensus/experimental/Filters.h>

using Interval = PacBio::Data::Interval;

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

std::vector<CoverageInterval> CoverageIntervals(const Interval& window,
                                                const std::vector<Interval>& input)
{
    std::vector<CoverageInterval> result;
    if (input.empty()) {
        result.emplace_back(CoverageInterval{window, 0});
        return result;
    }

    const auto winStart = window.Left();
    const auto winEnd = window.Right();
    auto pos = winStart;

    // leading uncovered interval
    CoverageInterval currentCoverageInterval{input.at(0), 1};
    if (pos < currentCoverageInterval.interval.Left()) {
        result.emplace_back(
            CoverageInterval{Interval{pos, currentCoverageInterval.interval.Left()}, 0});
        pos = currentCoverageInterval.interval.Left();
    }

    // find covered intervals & holes
    for (size_t i = 1; i < input.size(); ++i) {
        auto nextInterval = input.at(i);
        if (currentCoverageInterval.interval.Overlaps(nextInterval)) {  // overlap/adjacent
            currentCoverageInterval.interval = currentCoverageInterval.interval.Union(nextInterval);
            ++currentCoverageInterval.coverage;
            pos = currentCoverageInterval.interval.Right();
        } else {  // disjoint
            result.emplace_back(currentCoverageInterval);
            pos = currentCoverageInterval.interval.Right();
            result.emplace_back(CoverageInterval{Interval{pos, nextInterval.Left()}, 0});
            currentCoverageInterval = {nextInterval, 1};
            pos = currentCoverageInterval.interval.Right();
        }
    }

    // last covered interval
    result.emplace_back(currentCoverageInterval);
    pos = currentCoverageInterval.interval.Right();

    // last uncovered interval
    if (pos < winEnd) result.emplace_back(CoverageInterval{Interval{pos, winEnd}, 0});

    // clip first/last intervals to window bounds
    assert(!result.empty());
    result.front().interval = result.front().interval.Intersect(window);
    result.back().interval = result.back().interval.Intersect(window);
    return result;
}

std::vector<PacBio::Data::Interval> Holes(const PacBio::Data::Interval& windowInterval,
                                          const std::vector<PacBio::Data::Interval>& intervals)
{
    std::vector<PacBio::Data::Interval> result;

    size_t lastE = windowInterval.Left();
    for (const auto i : intervals) {
        if (i.Left() > lastE) result.push_back({lastE, i.Left()});
        lastE = i.Right();
    }
    if (lastE < windowInterval.Right()) result.push_back({lastE, windowInterval.Right()});

    return result;
}

std::vector<PacBio::Data::Interval> FancyIntervals(
    const PacBio::Data::Interval& windowInterval,
    const std::vector<PacBio::Data::Interval>& readIntervals, const size_t minCoverage)
{
    // transform readIntervals into extended intervals with minCoverage, and the
    // remainder as 'holes'
    auto intervals = KSpannedIntervals(windowInterval, std::move(readIntervals), minCoverage);
    auto holes = Holes(windowInterval, intervals);

    // concatenate, sort, & return intervals
    std::vector<PacBio::Data::Interval> result{std::move(intervals)};
    for (auto& hole : holes)
        result.emplace_back(std::move(hole));
    std::sort(result.begin(), result.end());
    return result;
}

std::vector<PacBio::Data::Interval> FancyIntervals(const PacBio::BAM::PbiRawData& index,
                                                   const ReferenceWindow& window,
                                                   const size_t minCoverage, const uint8_t minMapQV)
{
    const auto readIntervals = FilteredWindowIntervals(index, window, minMapQV);
    return FancyIntervals(window.interval, readIntervals, minCoverage);
}

std::vector<PacBio::Data::Interval> KSpannedIntervals(
    const PacBio::Data::Interval& windowInterval, std::vector<PacBio::Data::Interval> readIntervals,
    const size_t minCoverage, const size_t minLength)
{
    const auto k = minCoverage;
    assert(k >= 1);
    const auto winStart_ = windowInterval.Left();
    const auto winEnd_ = windowInterval.Right();

    // clip to window & translate to coordinates within window
    for (auto& interval : readIntervals) {
        interval = interval.Intersect(windowInterval);
        interval.Reset(interval.Left() - winStart_, interval.Right() - winStart_);
    }

    size_t winStart = 0;
    size_t winEnd = winEnd_ - winStart_;

    std::vector<size_t> positions(winEnd - winStart);
    std::iota(positions.begin(), positions.end(), 0);

    const auto coverage = ProjectIntoRange(readIntervals, windowInterval);
    assert(positions.size() == coverage.size());

    std::vector<PacBio::Data::Interval> intervalsFound;
    size_t x = std::string::npos;
    size_t y = 0;
    while (y < winEnd) {
        // Step 1: let x be the first pos >= y that is k-covered
        auto findFirstKCoveredPos = [&]() -> size_t {
            for (size_t i = y; i < positions.size(); ++i) {
                if (coverage.at(i) >= minCoverage) return positions.at(i);
            }
            return std::string::npos;
        };

        x = findFirstKCoveredPos();
        if (x == std::string::npos) break;

        // Step 2: extend the window [x, y) until [x, y) is no longer
        // k-spanned.  Do this by setting y to the k-th largest `end`
        // among reads covering x
        auto findEligibleEnds = [&]() -> std::vector<size_t> {
            std::vector<size_t> result;
            for (const auto& interval : readIntervals) {
                if (interval.Left() <= x) result.push_back(interval.Right());
            }
            return result;
        };

        auto eligibleEnds = findEligibleEnds();
        if (eligibleEnds.size() < minCoverage) break;
        std::sort(eligibleEnds.begin(), eligibleEnds.end());
        y = *(eligibleEnds.rbegin() + (minCoverage - 1));  // y = eligible[-k]

        // store extended interval
        intervalsFound.push_back({x, y});
    }

    // translate intervals back (respecting requesting minLength)
    std::vector<PacBio::Data::Interval> result;
    for (const auto& interval : intervalsFound) {
        if (interval.Length() >= minLength) {
            result.push_back({interval.Left() + winStart_, interval.Right() + winStart_});
        }
    }
    return result;
}

std::vector<PacBio::Data::Interval> FilteredIntervals(const PacBio::BAM::PbiRawData& index,
                                                      const PacBio::BAM::PbiFilter& filter)
{
    std::vector<PacBio::Data::Interval> readIntervals;

    const auto& mappedData = index.MappedData();
    for (size_t i = 0; i < index.NumReads(); ++i) {
        if (filter.Accepts(index, i)) {
            const auto start = mappedData.tStart_.at(i);
            const auto end = mappedData.tEnd_.at(i);
            readIntervals.emplace_back(start, end);
        }
    }

    return readIntervals;
}

std::vector<PacBio::Data::Interval> FilteredWindowIntervals(const PacBio::BAM::PbiRawData& index,
                                                            const ReferenceWindow& window,
                                                            const uint8_t minMapQV)
{
    const auto filter = MakeWindowFilter(window, minMapQV);
    auto readIntervals = FilteredIntervals(index, filter);
    std::sort(readIntervals.begin(), readIntervals.end());
    return readIntervals;
}

std::vector<size_t> ProjectIntoRange(const std::vector<PacBio::Data::Interval>& intervals,
                                     const PacBio::Data::Interval& windowInterval)
{
    std::vector<size_t> result(windowInterval.Length(), 0);

    const auto winStart = windowInterval.Left();
    const auto winEnd = windowInterval.Right();
    for (const auto& interval : intervals) {
        const size_t tStart = Clamp(interval.Left(), winStart, winEnd) - winStart;
        const size_t tEnd = Clamp(interval.Right(), winStart, winEnd) - winStart;
        for (size_t i = tStart; i < tEnd; ++i)
            ++result.at(i);
    }
    return result;
}

std::vector<PacBio::Data::Interval> SplitInterval(const PacBio::Data::Interval& source,
                                                  const size_t span)
{
    std::vector<Interval> result;

    size_t pos = source.Left();
    while (pos < source.Right()) {
        const auto left = pos;
        const auto right = pos + span;
        result.emplace_back(source.Intersect({left, right}));
        pos += span;
    }
    return result;
}

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
