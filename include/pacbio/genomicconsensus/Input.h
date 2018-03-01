// Author: Derek Barnett

#pragma once

#include <algorithm>
#include <cstdint>
#include <memory>
#include <tuple>
#include <vector>

#include <pbbam/BamRecord.h>
#include <pbbam/FastaSequenceQuery.h>
#include <pbbam/IndexedFastaReader.h>
#include <pbbam/PbiFilter.h>
#include <pbbam/PbiFilterQuery.h>

#include <pacbio/genomicconsensus/Intervals.h>
#include <pacbio/genomicconsensus/ReferenceWindow.h>
#include <pacbio/genomicconsensus/Settings.h>
#include <pacbio/genomicconsensus/Sorting.h>

namespace PacBio {
namespace GenomicConsensus {

class Input
{
public:
    explicit Input(const Settings& settings) : settings_{settings} {}

    Input() = delete;
    Input(const Input&) = delete;
    Input(Input&&) = default;
    Input& operator=(const Input&) = delete;
    Input& operator=(Input&&) = default;
    ~Input() = default;

public:
    ReferenceWindow EnlargedWindow(const ReferenceWindow& window) const;
    std::vector<PacBio::BAM::BamRecord> ReadsInWindow(const ReferenceWindow& window) const;
    std::string ReferenceInWindow(const ReferenceWindow& window) const;
    std::vector<ReferenceWindow> ReferenceWindows() const;

private:
    Settings settings_;
};

inline ReferenceWindow Input::EnlargedWindow(const ReferenceWindow& window) const
{
    PacBio::BAM::IndexedFastaReader fasta{settings_.referenceFilename};
    const PacBio::Data::Interval refInterval{
        0, static_cast<size_t>(fasta.SequenceLength(window.name))};

    const auto wStart = window.Start();
    const auto overhang = settings_.windowOverhang;

    const auto left = ((wStart < overhang) ? 0 : wStart - overhang);
    const auto right = window.End() + overhang;
    return ReferenceWindow{window.name, refInterval.Intersect({left, right})};
}

inline std::vector<PacBio::BAM::BamRecord> Input::ReadsInWindow(const ReferenceWindow& window) const
{
    using namespace PacBio::BAM;

    std::vector<BamRecord> result;
    std::vector<BamRecord> partialHits;
    result.reserve(settings_.maxCoverage);
    partialHits.reserve(settings_.maxCoverage * 2);

    const PbiFilter filter{
        {PbiReferenceEndFilter{static_cast<uint32_t>(window.Start()), Compare::GREATER_THAN},
         PbiReferenceStartFilter{static_cast<uint32_t>(window.End()), Compare::LESS_THAN},
         PbiMapQualityFilter{settings_.minMapQV, Compare::GREATER_THAN_EQUAL},
         PbiReferenceNameFilter{window.name}}};

    PbiFilterQuery query{filter, settings_.inputFilename};
    for (const auto& record : query) {
        const auto IsPoaCompatible = [&](const BamRecord& record) {
            const auto readLength = record.AlignedEnd() - record.AlignedStart();
            const auto refLength = record.ReferenceEnd() - record.ReferenceStart();
            const auto snr = record.SignalToNoise();
            return (readLength >= refLength * Settings::Defaults::ReadStumpinessThreshold) ||
                   (*std::min_element(snr.begin(), snr.end()) >= settings_.minHqRegionSnr) ||
                   (record.ReadAccuracy() >= settings_.minReadScore);
        };

        // quit if max coverage met
        if (result.size() == settings_.maxCoverage) break;

        // skip read if fails additional (non-PBI-backed) filters
        if (!IsPoaCompatible(record)) continue;

        // record spans window or exact hit
        const auto refStart = static_cast<uint32_t>(record.ReferenceStart());
        const auto refEnd = static_cast<uint32_t>(record.ReferenceEnd());
        if (refStart <= window.Start() && refEnd >= window.End())
            result.emplace_back(std::move(record));

        // read starts/ends within window
        else
            partialHits.emplace_back(std::move(record));
    }

    if (result.size() <= settings_.maxCoverage) {
        const auto partialHitLength = [&window](const BamRecord& r) {
            const auto refStart = static_cast<uint32_t>(r.ReferenceStart());
            const auto refEnd = static_cast<uint32_t>(r.ReferenceEnd());
            const auto winStart = window.Start();
            const auto winEnd = window.End();
            if (refStart > winStart)
                return winEnd - refStart;
            else
                return refEnd - winStart;
        };

        std::sort(partialHits.begin(), partialHits.end(),
                  [&partialHitLength](const BamRecord& lhs, const BamRecord& rhs) {
                      return partialHitLength(lhs) > partialHitLength(rhs);
                  });

        const auto maxNumHits = std::min(result.capacity() - result.size(), partialHits.size());
        for (size_t i = 0; i < maxNumHits; ++i)
            result.emplace_back(std::move(partialHits[i]));
    }

    Sorting::SortReadsInWindow(&result, window, settings_.sortStrategy);
    return result;
}

inline std::string Input::ReferenceInWindow(const ReferenceWindow& window) const
{
    PacBio::BAM::IndexedFastaReader reader{settings_.referenceFilename};
    return reader.Subsequence(window.name, static_cast<PacBio::BAM::Position>(window.Start()),
                              static_cast<PacBio::BAM::Position>(window.End()));
}

inline std::vector<ReferenceWindow> Input::ReferenceWindows() const
{
    using PacBio::Data::Interval;
    using PacBio::GenomicConsensus::SplitInterval;

    std::vector<ReferenceWindow> result;
    PacBio::BAM::FastaSequenceQuery query{settings_.referenceFilename};
    for (const auto& fasta : query) {
        const auto length = fasta.Bases().size();
        const auto source = Interval{0, length};
        const auto intervals = SplitInterval(source, settings_.windowSpan);

        const auto name = fasta.Name();
        for (const auto& interval : intervals)
            result.emplace_back(ReferenceWindow{name, interval});
    }
    return result;
}

}  // namespace GenomicConsensus
}  // namespace PacBio
