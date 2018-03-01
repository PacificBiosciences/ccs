// Author: Derek Barnett

#include <pacbio/genomicconsensus/experimental/Filters.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

void FilterAlignments(std::vector<BAM::BamRecord>* const reads, const float readStumpinessThreshold,
                      const float minHqRegionSnr, const float minReadScore)
{
    const auto IsIncompatible = [&](const PacBio::BAM::BamRecord& record) -> bool {
        const auto readLength = record.AlignedEnd() - record.AlignedStart();
        const auto refLength = record.ReferenceEnd() - record.ReferenceStart();
        const auto snr = record.SignalToNoise();
        return (readLength < refLength * readStumpinessThreshold) ||
               (*std::min_element(snr.begin(), snr.end()) < minHqRegionSnr) ||
               (record.ReadAccuracy() < minReadScore);
    };

    reads->erase(std::remove_if(reads->begin(), reads->end(), IsIncompatible), reads->end());
}

void FilterVariants(std::vector<Variant>* const variants, const size_t minCoverage,
                    const size_t minConfidence)
{
    const auto IsIncompatible = [&](const Variant& v) -> bool {
        return v.coverage.get() < minCoverage || v.confidence.get() < minConfidence;
    };

    variants->erase(std::remove_if(variants->begin(), variants->end(), IsIncompatible),
                    variants->end());
}

PacBio::BAM::PbiFilter MakeWindowFilter(const ReferenceWindow& window, const size_t refId,
                                        const uint8_t minMapQV)
{
    // TODO: barcode filter??

    using namespace PacBio::BAM;
    return PbiFilter{
        {PbiReferenceIdFilter{static_cast<int32_t>(refId)},
         PbiReferenceStartFilter{static_cast<uint32_t>(window.Start()),
                                 Compare::GREATER_THAN_EQUAL},
         PbiReferenceEndFilter{static_cast<uint32_t>(window.End()), Compare::LESS_THAN},
         PbiMapQualityFilter{minMapQV, Compare::GREATER_THAN_EQUAL}}};
}

PacBio::BAM::PbiFilter MakeWindowFilter(const ReferenceWindow& window, const uint8_t minMapQV)
{
    // TODO: barcode filter??

    using namespace PacBio::BAM;
    return PbiFilter{
        {PbiReferenceEndFilter{static_cast<uint32_t>(window.Start()), Compare::GREATER_THAN},
         PbiReferenceStartFilter{static_cast<uint32_t>(window.End()), Compare::LESS_THAN},
         PbiMapQualityFilter{minMapQV, Compare::GREATER_THAN_EQUAL},
         PbiReferenceNameFilter{window.name}}};
}

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
