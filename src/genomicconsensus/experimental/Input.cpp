// Author: Derek Barnett

#include <pacbio/genomicconsensus/experimental/Input.h>

#include <pbbam/FastaSequenceQuery.h>
#include <pbbam/PbiFilterQuery.h>

#include <pacbio/genomicconsensus/experimental/Filters.h>
#include <pacbio/genomicconsensus/experimental/Intervals.h>
#include <pacbio/genomicconsensus/experimental/Sorting.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

Input::Input(const Settings& settings) : settings_{settings}, fasta_{settings.referenceFilename} {}

std::vector<PacBio::BAM::BamRecord> Input::ReadsInWindow(const ReferenceWindow& window) const
{
    using namespace PacBio::BAM;

    std::vector<BamRecord> result;
    std::vector<BamRecord> partialHits;
    result.reserve(settings_.maxCoverage);
    partialHits.reserve(settings_.maxCoverage * 2);

    const auto filter = MakeWindowFilter(window, settings_);
    PbiFilterQuery query{filter, settings_.inputFilename};
    for (const auto& record : query) {

        // TODO (DB): combine with similar lambda used in FilterAlignments (Filters.cpp),
        //            if possible?
        const auto IsPoaCompatible = [&](const BamRecord& record) {
            const auto readLength = record.AlignedEnd() - record.AlignedStart();
            const auto refLength = record.ReferenceEnd() - record.ReferenceStart();
            const auto snr = record.SignalToNoise();
            return (readLength >= refLength * settings_.readStumpinessThreshold) &&
                   (*std::min_element(snr.begin(), snr.end()) >= settings_.minHqRegionSnr) &&
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
                  [&window, &partialHitLength](const BamRecord& lhs, const BamRecord& rhs) {
                      return partialHitLength(lhs) > partialHitLength(rhs);
                  });

        const auto maxNumHits = std::min(result.capacity() - result.size(), partialHits.size());
        for (size_t i = 0; i < maxNumHits; ++i)
            result.emplace_back(std::move(partialHits[i]));
    }

    SortReadsInWindow(&result, window, settings_.sortStrategy);
    return result;
}

std::string Input::ReferenceInWindow(const ReferenceWindow& window) const
{
    return fasta_.Subsequence(window.name, static_cast<PacBio::BAM::Position>(window.Start()),
                              static_cast<PacBio::BAM::Position>(window.End()));
}

std::vector<std::string> Input::ReferenceNames() const
{
    std::vector<std::string> result;
    PacBio::BAM::FastaSequenceQuery query{settings_.referenceFilename};
    for (const auto& fasta : query)
        result.emplace_back(fasta.Name());
    return result;
}

std::vector<ReferenceWindow> Input::ReferenceWindows(bool splitWindows) const
{
    using PacBio::Data::Interval;
    using PacBio::GenomicConsensus::experimental::SplitInterval;

    std::vector<ReferenceWindow> result;
    PacBio::BAM::FastaSequenceQuery query{settings_.referenceFilename};
    for (const auto& fasta : query) {
        const auto name = fasta.Name();
        const auto length = fasta.Bases().size();
        const auto source = Interval{0, length};

        std::vector<PacBio::Data::Interval> intervals;
        if (splitWindows)
            intervals = SplitInterval(source, settings_.windowSpan);
        else
            intervals.push_back(source);

        for (const auto& interval : intervals)
            result.emplace_back(ReferenceWindow{name, interval});
    }
    return result;
}

size_t Input::SequenceLength(const std::string& refName) const
{
    return fasta_.SequenceLength(refName);
}

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
