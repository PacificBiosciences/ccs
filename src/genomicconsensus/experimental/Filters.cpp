// Copyright (c) 2017, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

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
