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

#pragma once

#include <vector>

#include <pbbam/BamRecord.h>
#include <pbbam/PbiFilter.h>

#include <pacbio/genomicconsensus/experimental/Settings.h>
#include <pacbio/genomicconsensus/experimental/Variant.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

///
/// \brief FilterAlignments
///
/// Filter in-place.
///
/// \param reads
/// \param readStumpinessThreshold
/// \param minHqRegionSnr
/// \param minReadScore
///
void FilterAlignments(std::vector<PacBio::BAM::BamRecord>* const reads,
                      const float readStumpinessThreshold, const float minHqRegionSnr,
                      const float minReadScore);

///
/// \brief FilterAlignments
///
/// Overloaded for Settings. Filter in-place.
///
/// \param reads
/// \param settings
///
inline void FilterAlignments(std::vector<PacBio::BAM::BamRecord>* const reads,
                             const Settings& settings)
{
    return FilterAlignments(reads, settings.readStumpinessThreshold, settings.minHqRegionSnr,
                            settings.minReadScore);
}

///
/// \brief FilteredAlignments
///
/// Return filtered copy.
///
/// \param reads
/// \param readStumpinessThreshold
/// \param minHqRegionSnr
/// \param minReadScore
/// \return
///
inline std::vector<PacBio::BAM::BamRecord> FilteredAlignments(
    const std::vector<PacBio::BAM::BamRecord>& reads, const float readStumpinessThreshold,
    const float minHqRegionSnr, const float minReadScore)
{
    auto v = reads;
    FilterAlignments(&v, readStumpinessThreshold, minHqRegionSnr, minReadScore);
    return v;
}

///
/// \brief FilteredAlignments
///
/// Overloaded for Settings. Return filtered copy.
///
/// \param reads
/// \param settings
/// \return
///
inline std::vector<PacBio::BAM::BamRecord> FilteredAlignments(
    const std::vector<PacBio::BAM::BamRecord>& reads, const Settings& settings)
{
    auto v = reads;
    FilterAlignments(&v, settings.readStumpinessThreshold, settings.minHqRegionSnr,
                     settings.minReadScore);
    return v;
}

///
/// \brief FilterVariants
///
/// Filter in-place.
///
/// \param variants
/// \param minCoverage
/// \param minConfidence
///
void FilterVariants(std::vector<Variant>* const variants, const size_t minCoverage,
                    const size_t minConfidence);

///
/// \brief FilterVariants
///
/// Overloaded for Settings. Filter in-place.
///
/// \param variants
/// \param settings
///
inline void FilterVariants(std::vector<Variant>* variants, const Settings& settings)
{
    return FilterVariants(variants, settings.minCoverage, settings.minConfidence);
}

///
/// \brief FilteredVariants
///
/// Return filtered copy.
///
/// \param variants
/// \param minCoverage
/// \param minConfidence
/// \return
///
inline std::vector<Variant> FilteredVariants(const std::vector<Variant>& variants,
                                             const size_t minCoverage, const size_t minConfidence)
{
    auto v = variants;
    FilterVariants(&v, minCoverage, minConfidence);
    return v;
}

///
/// \brief FilteredVariants
///
/// Overloaded for Settings. Return filtered copy.
///
/// \param variants
/// \param settings
/// \return
///
inline std::vector<Variant> FilteredVariants(const std::vector<Variant>& variants,
                                             const Settings& settings)
{
    return FilteredVariants(variants, settings.minCoverage, settings.minConfidence);
}

///
/// \brief MakeWindowFilter
///
/// Makes PbiFilter on window, with a minimum mapQV. Filtering using refId,
/// if available, is more efficient than using refName (window.name)
///
/// \param window
/// \param refId
/// \param minMapQV
/// \return
///
PacBio::BAM::PbiFilter MakeWindowFilter(const ReferenceWindow& window, const size_t refId,
                                        const uint8_t minMapQV);

///
/// \brief MakeWindowFilter
///
/// Overloaded for settings.
///
/// Makes PbiFilter on window, with a minimum mapQV. Filtering using refId,
/// if available, is more efficient than using refName (window.name)
///
/// \param window
/// \param refId
/// \param minMapQV
/// \return
///
inline PacBio::BAM::PbiFilter MakeWindowFilter(const ReferenceWindow& window, const size_t refId,
                                               const Settings& settings)
{
    return MakeWindowFilter(window, refId, settings.minMapQV);
}

///
/// \brief MakeWindowFilter
///
/// Makes PbiFilter on window, with a minimum mapQV.
///
/// \param window
/// \param minMapQV
/// \return
///
PacBio::BAM::PbiFilter MakeWindowFilter(const ReferenceWindow& window, const uint8_t minMapQV);

///
/// \brief MakeWindowFilter
///
/// Overloaded for Settings. Makes PbiFilter on window, with a minimum mapQV.
///
/// \param window
/// \param settings
/// \return
///
inline PacBio::BAM::PbiFilter MakeWindowFilter(const ReferenceWindow& window,
                                               const Settings& settings)
{
    return MakeWindowFilter(window, settings.minMapQV);
}

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
