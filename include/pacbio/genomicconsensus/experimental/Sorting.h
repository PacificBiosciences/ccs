// Author: Derek Barnett

#pragma once

#include <vector>

#include <pbbam/BamRecord.h>

#include <pacbio/genomicconsensus/experimental/SortingStrategy.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

struct ReferenceWindow;

///
/// \brief SortReadsInWindow
///
/// Sort in-place.
///
/// \param reads
/// \param window
/// \param strategy
///
void SortReadsInWindow(std::vector<PacBio::BAM::BamRecord>* const reads,
                       const ReferenceWindow& window, SortingStrategy strategy);

///
/// \brief SortedReadsInWindow
///
/// Return sorted copy.
///
/// \param reads
/// \param window
/// \param strategy
/// \return
///
inline std::vector<PacBio::BAM::BamRecord> SortedReadsInWindow(
    const std::vector<PacBio::BAM::BamRecord>& reads, const ReferenceWindow& window,
    SortingStrategy strategy)
{
    std::vector<PacBio::BAM::BamRecord> result = reads;
    SortReadsInWindow(&result, window, strategy);
    return result;
}

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
