// Author: Derek Barnett

#pragma once

#include <memory>
#include <stdexcept>
#include <tuple>
#include <vector>

#include <pbbam/BamRecord.h>

#include <pacbio/data/Read.h>
#include <pacbio/genomicconsensus/ReferenceWindow.h>
#include <pacbio/genomicconsensus/SortingStrategy.h>

namespace PacBio {
namespace GenomicConsensus {

struct Sorting
{
public:
    static void SortReadsInWindow(std::vector<PacBio::BAM::BamRecord>* reads,
                                  const ReferenceWindow& window, SortingStrategy strategy);

    static std::vector<PacBio::BAM::BamRecord> SortReadsInWindow(
        const std::vector<PacBio::BAM::BamRecord>& reads, const ReferenceWindow& window,
        SortingStrategy strategy);
};

inline void Sorting::SortReadsInWindow(std::vector<PacBio::BAM::BamRecord>* reads,
                                       const ReferenceWindow& window,
                                       const SortingStrategy strategy)
{
    using BamRecord = PacBio::BAM::BamRecord;

    auto lengthInWindow = [&window](const BamRecord& read) {
        const auto start = std::max(static_cast<size_t>(read.ReferenceStart()), window.Start());
        const auto end = std::min(static_cast<size_t>(read.ReferenceEnd()), window.End());
        return end - start;
    };

    switch (strategy) {
        case SortingStrategy::LONGEST_AND_STRAND_BALANCED: {
            // lexsort by read start, then end
            std::stable_sort(reads->begin(), reads->end(),
                             [](const BamRecord& lhs, const BamRecord& rhs) {
                                 const auto lhsStart = lhs.ReferenceStart();
                                 const auto rhsStart = rhs.ReferenceStart();
                                 const auto lhsEnd = lhs.ReferenceEnd();
                                 const auto rhsEnd = rhs.ReferenceEnd();
                                 return std::tie(lhsStart, lhsEnd) < std::tie(rhsStart, rhsEnd);
                             });

            // sort (descending) by read length
            std::stable_sort(reads->begin(), reads->end(),
                             [lengthInWindow](const BamRecord& lhs, const BamRecord& rhs) {
                                 const auto lhsLength = lengthInWindow(lhs);
                                 const auto rhsLength = lengthInWindow(rhs);
                                 return lhsLength > rhsLength;
                             });
            break;
        }
        case SortingStrategy::LONGEST: {
            std::stable_sort(reads->begin(), reads->end(),
                             [lengthInWindow](const PacBio::BAM::BamRecord& lhs,
                                              const PacBio::BAM::BamRecord& rhs) {
                                 const auto lhsLength = lengthInWindow(lhs);
                                 const auto rhsLength = lengthInWindow(rhs);
                                 return rhsLength < lhsLength;
                             });
            break;
        }
        case SortingStrategy::SPANNING: {
            reads->erase(
                std::remove_if(reads->begin(), reads->end(),
                               [&window, lengthInWindow](const PacBio::BAM::BamRecord& record) {
                                   const auto length = lengthInWindow(record);
                                   return length != window.Length();
                               }),
                reads->end());
            break;
        }
        case SortingStrategy::FILE_ORDER: {
            // no sorting necessary
            break;
        }
        default:
            throw std::runtime_error("unexpected sorting strategy");
    }
}

inline std::vector<PacBio::BAM::BamRecord> Sorting::SortReadsInWindow(
    const std::vector<PacBio::BAM::BamRecord>& reads, const ReferenceWindow& window,
    const SortingStrategy strategy)
{
    std::vector<PacBio::BAM::BamRecord> result = reads;
    SortReadsInWindow(&result, window, strategy);
    return result;
}

}  // namespace GenomicConsensus
}  // namespace PacBio
