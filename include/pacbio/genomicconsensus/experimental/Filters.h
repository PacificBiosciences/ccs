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
