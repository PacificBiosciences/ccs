// Author: Derek Barnett

#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include <pbbam/BamRecord.h>

#include <pacbio/UnanimityConfig.h>

#include <pacbio/align/AffineAlignment.h>
#include <pacbio/align/AlignConfig.h>
#include <pacbio/align/PairwiseAlignment.h>
#include <pacbio/data/Interval.h>
#include <pacbio/data/internal/BaseEncoding.h>
#include <pacbio/denovo/PoaConsensus.h>
#include <pacbio/genomicconsensus/experimental/ReferenceWindow.h>
#include <pacbio/genomicconsensus/experimental/Settings.h>
#include <pacbio/genomicconsensus/experimental/Variant.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

// NOTE: This are structs/methods common between Arrow & POA models. Upcoming PR
//       will switch Arrow over to using these. Likely candidate for refactoring
//       and/or renaming but we're building & iterating up to that point.

template <typename T>
static inline bool FoundCode(const T& validCodes, const char code)
{
    for (const char c : validCodes) {
        if (c == code) return true;
    }
    return false;
}

static inline void ClipReadsToWindow(std::vector<PacBio::BAM::BamRecord>* const reads,
                                     const ReferenceWindow& window)
{
    const auto winStart = window.Start();
    const auto winEnd = window.End();
    for (auto& read : *reads)
        read.Clip(PacBio::BAM::ClipType::CLIP_TO_REFERENCE, winStart, winEnd);
}

static inline void FilterAlignments(std::vector<PacBio::BAM::BamRecord>* const reads,
                                    const Settings& settings)
{
    const auto IsPoaIncompatible = [&](const PacBio::BAM::BamRecord& record) -> bool {
        const auto readLength = record.AlignedEnd() - record.AlignedStart();
        const auto refLength = record.ReferenceEnd() - record.ReferenceStart();
        const auto snr = record.SignalToNoise();
        return (readLength < refLength * settings.readStumpinessThreshold) ||
               (*std::min_element(snr.begin(), snr.end()) < settings.minHqRegionSnr) ||
               (record.ReadAccuracy() < settings.minReadScore);
    };

    reads->erase(std::remove_if(reads->begin(), reads->end(), IsPoaIncompatible), reads->end());
}

static inline std::vector<std::string> FilteredForwardSequences(
    const std::vector<PacBio::BAM::BamRecord>& reads, const ReferenceWindow& window)
{
    using BamRecord = PacBio::BAM::BamRecord;
    using Orientation = PacBio::BAM::Orientation;
    using Strand = PacBio::BAM::Strand;

    const auto spansReferenceRange = [](const BamRecord& read, const ReferenceWindow& window) {
        const auto tStart = static_cast<size_t>(read.ReferenceStart());
        const auto tEnd = static_cast<size_t>(read.ReferenceEnd());
        assert(window.Start() <= window.End());
        return (tStart <= window.Start() && tEnd >= window.End());
    };

    std::vector<std::string> result;
    for (const auto& read : reads) {
        if (read.AlignedStrand() == Strand::FORWARD && spansReferenceRange(read, window))
            result.push_back(read.Sequence(Orientation::NATIVE, false));  // excise soft clips?
    }
    return result;
}

static std::vector<Variant> FilterVariants(const std::vector<Variant>& variants,
                                           const Settings& settings)
{
    std::vector<Variant> result;
    for (const auto& v : variants) {
        if (v.coverage.get() >= settings.minCoverage &&
            v.confidence.get() >= settings.minConfidence) {
            result.emplace_back(v);
        }
    }
    return result;
}

static void AnnotateVariants(std::vector<Variant>* const variants,
                             const std::vector<PacBio::BAM::BamRecord>& reads)
{
    for (auto& v : *variants) {
        std::string annotation;
        for (const auto& read : reads) {
            if (!annotation.empty()) {
                annotation.push_back(',');
                annotation.push_back(' ');
            }
            annotation.append(read.FullName());
        }
        v.Annotate("rows", annotation);
    }
}

static size_t Median(std::vector<size_t> v)
{
    const auto mid = v.size() / 2;
    auto nth = v.begin() + mid;
    std::nth_element(v.begin(), nth, v.end());
    if (v.size() % 2 == 0) {
        // avg of middle two
        const auto first = *nth;
        std::nth_element(v.begin(), --nth, v.end());
        const auto second = *nth;
        return (first + second) / 2;
    }
    return *nth;
}

static std::unique_ptr<const PacBio::Poa::PoaConsensus> MakePoaConsensus(
    std::vector<std::string>&& fwdSequences, const Settings& settings)
{
    using AlignMode = PacBio::Align::AlignMode;
    using PoaConsensus = PacBio::Poa::PoaConsensus;

    std::vector<size_t> seqLengths;
    seqLengths.reserve(fwdSequences.size());
    for (const auto& seq : fwdSequences)
        seqLengths.push_back(seq.size());

    const auto median = Median(std::move(seqLengths));
    std::vector<std::string> ordSeqs = std::move(fwdSequences);
    std::sort(ordSeqs.begin(), ordSeqs.end(),
              [median](const std::string& lhs, const std::string& rhs) {
                  const auto lhsSize = static_cast<int>(lhs.size());
                  const auto rhsSize = static_cast<int>(rhs.size());
                  const auto intMedian = static_cast<int>(median);
                  const auto diffLeft = std::abs(lhsSize - intMedian);
                  const auto diffRight = std::abs(rhsSize - intMedian);
                  return diffLeft < diffRight;
              });
    ordSeqs.resize(std::distance(ordSeqs.begin(), ordSeqs.begin() + settings.maxPoaCoverage));

    const auto poaConfig = PacBio::Poa::DefaultPoaConfig(AlignMode::GLOBAL);
    const auto cov = ordSeqs.size();
    const auto minCov = (cov < 5 ? 1 : ((cov + 1) / 2 - 1));
    return std::unique_ptr<const PacBio::Poa::PoaConsensus>{
        PoaConsensus::FindConsensus(ordSeqs, poaConfig, minCov)};
}

static std::vector<PacBio::Data::Interval> TranscriptIntervals(const std::string& transcript)
{
    using Interval = PacBio::Data::Interval;

    std::vector<Interval> result;
    const auto tSize = transcript.size();
    if (tSize == 0) return result;
    assert(!transcript.empty());

    char previousChar = transcript.at(0);
    size_t currentRunStart = 0;
    size_t currentRunLength = 1;

    for (size_t i = 1; i < tSize; ++i) {
        const auto currentChar = transcript.at(i);
        if (currentChar == previousChar)
            ++currentRunLength;
        else {
            result.emplace_back(currentRunStart, currentRunStart + currentRunLength);
            currentRunStart = i;
            currentRunLength = 1;
        }
        previousChar = currentChar;
    }
    // push last interval & return
    result.emplace_back(currentRunStart, currentRunStart + currentRunLength);
    return result;
}

static UNANIMITY_CONSTEXPR const char* LookupIUPAC(const char c)
{
    constexpr const std::array<const char*, 16> table{{"-", "A", "C", "AC", "G", "AG", "CG", "ACG",
                                                       "T", "AT", "CT", "ACT", "GT", "AGT", "CGT",
                                                       "ACGT"}};
    const auto ncbi4na = PacBio::Data::detail::NCBI4na::FromASCII(c);
    return table.at(ncbi4na.Data());
}

struct SplitupIUPACResult
{
    std::string readSeq1;
    boost::optional<std::string> readSeq2 = boost::none;
    boost::optional<double> freq1 = boost::none;
    boost::optional<double> freq2 = boost::none;
};

static SplitupIUPACResult SplitupIUPAC(const std::string& css)
{
    std::string listSeq1;
    std::string listSeq2;
    listSeq1.reserve(css.size());
    listSeq2.reserve(css.size());
    for (const auto c : css) {
        const std::string value{LookupIUPAC(c)};
        listSeq1.push_back(value.front());
        listSeq2.push_back(value.back());
    }

    SplitupIUPACResult result;
    if (listSeq1 == listSeq2) {
        // haploid
        result.readSeq1 = std::move(listSeq1);
    } else {
        // diploid
        result.readSeq1 = std::move(listSeq1);
        result.readSeq2 = std::move(listSeq2);
        result.freq1 = 0.5;
        result.freq2 = 0.5;
    }
    return result;
}

static std::vector<Variant> VariantsFromAlignment(
    const std::unique_ptr<PacBio::Align::PairwiseAlignment>& alignment,
    const ReferenceWindow& window, const boost::optional<std::vector<uint8_t>>& cssQvInWindow,
    const std::vector<uint8_t>& siteCoverage,
    const boost::optional<std::vector<uint8_t>>& effectiveSiteCoverage)
{
    std::vector<Variant> variants;

    const auto refId = window.name;
    const auto refStart = window.Start();
    size_t refPos = refStart;
    size_t cssPos = 0;
    char refPrev = 'N';
    char cssPrev = 'N';

    constexpr const std::array<char, 5> validCodes{{'R', 'I', 'D', 'M', 'N'}};

    // We don't call variants where either the reference or css is 'N'
    const auto& target = alignment->Target();
    const auto& query = alignment->Query();
    auto transcript = alignment->Transcript();
    for (size_t i = 0; i < transcript.size(); ++i) {
        const char t = target.at(i);
        const char q = query.at(i);
        if (t == 'N' || q == 'N') transcript[i] = 'N';
    }

    boost::optional<Variant> variant;
    for (const auto& interval : TranscriptIntervals(transcript)) {

        const auto pos = interval.Left();
        const auto code = transcript.at(pos);
        if (!FoundCode(validCodes, code)) throw std::runtime_error{"invalid code"};

        const auto length = interval.Length();
        auto ref = target.substr(pos, length);
        auto css = query.substr(pos, length);

        const auto NotGap = [](const char c) { return c != '-'; };
        const auto refLen = std::count_if(ref.begin(), ref.end(), NotGap);
        const auto cssLen = std::count_if(ref.begin(), ref.end(), NotGap);

        variant = boost::none;

        switch (code) {
            case 'M':  // fall-through
            case 'N':  // .
                break;

            case 'R': {
                assert(css.size() == ref.size());
                const auto splitupResult = SplitupIUPAC(css);
                css = splitupResult.readSeq1;
                variant = Variant{refId,   refPos, refPos + css.size(), ref, splitupResult.readSeq1,
                                  refPrev, cssPrev};
                variant->readSeq2 = splitupResult.readSeq2;
                variant->frequency1 = splitupResult.freq1;
                variant->frequency2 = splitupResult.freq2;
                break;
            }

            case 'I': {
                const auto splitupResult = SplitupIUPAC(css);
                css = splitupResult.readSeq1;
                variant =
                    Variant{refId, refPos, refPos, "", splitupResult.readSeq1, refPrev, cssPrev};
                variant->readSeq2 = splitupResult.readSeq2;
                variant->frequency1 = splitupResult.freq1;
                variant->frequency2 = splitupResult.freq2;
                break;
            }

            case 'D': {
                variant = Variant{refId, refPos, refPos + ref.size(), ref, "", refPrev, cssPrev};
                break;
            }

            default:
                throw std::runtime_error("invalid code");
        }

        if (variant) {
            // HACK ALERT: variants at the first and last position
            // are not handled correctly

            if (!siteCoverage.empty()) {
                const auto i = std::min(refPos - window.Start(), siteCoverage.size() - 1);
                variant->coverage = siteCoverage.at(i);
            }
            if (effectiveSiteCoverage && !effectiveSiteCoverage->empty()) {
                const auto i = std::min(refPos - window.Start(), effectiveSiteCoverage->size() - 1);
                variant->Annotate("effectiveSiteoverage",
                                  std::to_string(effectiveSiteCoverage->at(i)));
            }
            if (cssQvInWindow && !cssQvInWindow->empty()) {
                const auto i = std::min(cssPos, cssQvInWindow->size() - 1);
                variant->confidence = cssQvInWindow->at(i);
            }
            variants.push_back(variant.get());
        }

        // update counters
        refPos += refLen;
        cssPos += cssLen;

        const auto IsGap = [](const char c) { return c == '-'; };
        ref.erase(std::remove_if(ref.begin(), ref.end(), IsGap), ref.end());
        css.erase(std::remove_if(css.begin(), css.end(), IsGap), css.end());

        refPrev = (ref.empty() ? refPrev : ref.back());
        cssPrev = (css.empty() ? cssPrev : css.back());
    }
    return variants;
}

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
