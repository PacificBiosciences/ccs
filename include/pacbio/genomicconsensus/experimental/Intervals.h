// Author: Derek Barnett

#pragma once

#include <cstddef>
#include <vector>

#include <pacbio/data/Interval.h>

#include <pacbio/genomicconsensus/experimental/Settings.h>

namespace PacBio {
namespace BAM {

class PbiFilter;
class PbiRawData;

}  // namespace BAM
}  // namespace PacBio

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

struct ReferenceWindow;

///
/// \brief Clamp
/// \param pos
/// \param min
/// \param max
/// \return
///
inline size_t Clamp(const size_t pos, const size_t min, const size_t max)
{
    // TODO (DB): Probably won't leave this here. Maybe a more general math
    // utils or somthing. But this used by Intervals.cpp, so it's fine for now.

    if (pos < min) return min;
    if (pos > max) return max;
    return pos;
}

///
/// \brief The CoverageInterval struct
///
struct CoverageInterval
{
    PacBio::Data::Interval interval;
    size_t coverage;
};

///
/// \brief CoverageIntervals
/// \param window
/// \param input
/// \return
///
std::vector<CoverageInterval> CoverageIntervals(const PacBio::Data::Interval& window,
                                                const std::vector<PacBio::Data::Interval>& input);

///
/// \brief FancyIntervals
///
///
///
/// \param window
/// \param readIntervals
/// \param minCoverage
/// \return
///
std::vector<PacBio::Data::Interval> FancyIntervals(
    const PacBio::Data::Interval& windowInterval,
    const std::vector<PacBio::Data::Interval>& readIntervals, const size_t minCoverage);

///
/// \brief FancyIntervals
///
/// Finds a maximal set of maximal disjoint intervals within \p window such that
/// each interval is spanned by at least minCoverage reads.
///
/// Note that this is a greedy search procedure and may not always return the
/// optimal solution, in some sense.  However it will always return the optimal
/// solutions in the most common cases.
///
/// Fills in the remaining gaps, and adds them to output.
///
std::vector<PacBio::Data::Interval> FancyIntervals(const PacBio::BAM::PbiRawData& index,
                                                   const ReferenceWindow& window,
                                                   const size_t minCoverage,
                                                   const uint8_t minMapQV);

///
/// \brief FancyIntervals
///
/// Override for Settings. Finds a maximal set of maximal disjoint intervals
/// within \p window such that each interval is spanned by at least minCoverage
/// reads.
///
/// \param index
/// \param window
/// \param settings
/// \return
///
inline std::vector<PacBio::Data::Interval> FancyIntervals(const PacBio::BAM::PbiRawData& index,
                                                          const ReferenceWindow& window,
                                                          const Settings& settings)
{
    return FancyIntervals(index, window, settings.minCoverage, settings.minMapQV);
}

/// \brief FilteredIntervals
///
/// \param index  PBI index data
/// \param filter criteria to filter on
/// \return
///
std::vector<PacBio::Data::Interval> FilteredIntervals(const PacBio::BAM::PbiRawData& index,
                                                      const PacBio::BAM::PbiFilter& filter);

///
/// \brief FilteredWindowIntervals
///
/// Return sorted read intervals within window, satisfying minMapQV. Equivalent to:
///
/// \param index
/// \param window
/// \param minMapQV
/// \return
///
std::vector<PacBio::Data::Interval> FilteredWindowIntervals(const PacBio::BAM::PbiRawData& index,
                                                            const ReferenceWindow& window,
                                                            const uint8_t minMapQV);

/////
///// \brief Holes
/////
///// Given a window and a set of disjoint subintervals, return the "holes", which
///// are the intervals of the refWindow not covered by the given subintervals.
/////
///// \param window
///// \param intervals
///// \return
/////
std::vector<PacBio::Data::Interval> Holes(const Data::Interval& windowInterval,
                                          const std::vector<PacBio::Data::Interval>& intervals);

/////
///// \brief KSpannedIntervals
/////
///// Find a maximal set of maximal disjoint intervals within \p window such that
///// each interval is spanned by at least \p minCoverage reads.
/////
///// Note that this is a greedy search procedure and may not always return the
///// optimal solution, in some sense. However it will always return the optimal
///// solutions in the most common cases.
/////
///// \param windowInterval     window under consideration
///// \param readIntervals      read intervals, sorted in ascending order.
///// \param minCoverage        the number of reads that must span intervals to be returned
/////
///// \return sorted list of intervals
/////
std::vector<PacBio::Data::Interval> KSpannedIntervals(
    const PacBio::Data::Interval& windowInterval, std::vector<PacBio::Data::Interval> readIntervals,
    const size_t minCoverage, const size_t minLength = 0);

///
/// \brief ProjectIntoRange
///
/// Find coverage in the \p window implied by \p  intervals.
///
/// \param intervals
/// \param windowInterval
/// \return
///
std::vector<size_t> ProjectIntoRange(const std::vector<PacBio::Data::Interval>& intervals,
                                     const Data::Interval& windowInterval);

///
/// \brief SplitInterval
/// \param source
/// \param span
/// \return
///
std::vector<PacBio::Data::Interval> SplitInterval(const PacBio::Data::Interval& source,
                                                  const size_t span);

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
