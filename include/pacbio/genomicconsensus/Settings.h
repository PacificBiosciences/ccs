// Author: Derek Barnett

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include <pacbio/genomicconsensus/SortingStrategy.h>

namespace PacBio {
namespace GenomicConsensus {

struct Settings
{
    struct Defaults
    {
        // CLI-exposed
        static constexpr size_t WindowSpan = 500;
        static constexpr size_t WindowOverhang = 5;
        static constexpr uint8_t MinMapQV = 10;
        static constexpr size_t MinCoverage = 5;
        static constexpr size_t MaxCoverage = 100;
        static constexpr float MinReadScore = 0.65;
        static constexpr float MinHqRegionSnr = 3.75;
        static constexpr size_t MinPoaCoverage = 3;
        static constexpr size_t MaxPoaCoverage = 11;
        static constexpr SortingStrategy Strategy = SortingStrategy::LONGEST_AND_STRAND_BALANCED;
        static constexpr float ReadStumpinessThreshold = 0.1;
        static constexpr bool AnnotateGFF = false;
        static constexpr bool ReportEffectiveCoverage = false;
        static constexpr bool PolishDiploid = true;
        static constexpr size_t MinConfidence = 40;
        static constexpr double MinZScore = -3.4;
        static constexpr bool ComputeConfidence = true;
        static constexpr double MinAccuracy = 0.82;
        static constexpr size_t MaxIterations = 40;
        static constexpr size_t MutationSeparation = 10;
        static constexpr size_t MutationNeighborhood = 20;
        static constexpr size_t MaskRadius = 0;
        static constexpr double MaskErrorRate = 0.0;
        static constexpr bool Diploid = false;
    };

    std::string inputFilename;
    std::string referenceFilename;
    std::string fastaFilename;
    std::string fastqFilename;
    std::string gffFilename;
    std::string vcfFilename;

    size_t windowSpan = Defaults::WindowSpan;
    size_t windowOverhang = Defaults::WindowOverhang;
    uint8_t minMapQV = Defaults::MinMapQV;
    size_t minCoverage = Defaults::MinCoverage;
    size_t maxCoverage = Defaults::MaxCoverage;
    float minReadScore = Defaults::MinReadScore;
    float minHqRegionSnr = Defaults::MinHqRegionSnr;
    SortingStrategy sortStrategy = Defaults::Strategy;
    float readStumpinessThreshold = Defaults::ReadStumpinessThreshold;
    size_t minPoaCoverage = Defaults::MinPoaCoverage;
    size_t maxPoaCoverage = Defaults::MaxPoaCoverage;
    size_t minConfidence = Defaults::MinConfidence;
    bool annotateGFF = Defaults::AnnotateGFF;
    bool reportEffectiveCoverage = Defaults::ReportEffectiveCoverage;
    bool polishDiploid = Defaults::PolishDiploid;
    double minZScore = Defaults::MinZScore;
    bool computeConfidence = Defaults::ComputeConfidence;
    double minAccuracy = Defaults::MinAccuracy;
    size_t maxIterations = Defaults::MaxIterations;
    size_t mutationSeparation = Defaults::MutationSeparation;
    size_t mutationNeighborhood = Defaults::MutationNeighborhood;
    size_t maskRadius = Defaults::MaskRadius;
    double maskErrorRate = Defaults::MaskErrorRate;
    bool diploid = Defaults::Diploid;
};

}  // namespace GenomicConsensus
}  // namespace PacBio
