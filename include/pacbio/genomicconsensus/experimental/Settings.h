// Author: Derek Barnett

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include <pacbio/genomicconsensus/experimental/ConsensusMode.h>
#include <pacbio/genomicconsensus/experimental/NoCallStyle.h>
#include <pacbio/genomicconsensus/experimental/ReferenceWindow.h>
#include <pacbio/genomicconsensus/experimental/SortingStrategy.h>
#include <pbcopper/cli/Results.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

///
/// \brief The Settings struct
///
struct Settings
{
    struct Defaults
    {
        // parallelism
        static constexpr const size_t NumThreads = 1;

        // output filtering
        static constexpr const size_t MinConfidence = 40;
        static constexpr const size_t MinCoverage = 5;
        static constexpr const NoCallStyle NoCall = NoCallStyle::NO_CALL;

        // read selection/filtering
        static constexpr const double MinAccuracy = 0.82;
        static constexpr const size_t MaxCoverage = 100;
        static constexpr const uint8_t MinMapQV = 10;
        static constexpr const float MinReadScore = 0.65;
        static constexpr const float MinHqRegionSnr = 3.75;
        static constexpr const double MinZScore = -3.4;

        // algorithm and parameters
        static constexpr const ConsensusMode Mode = ConsensusMode::ARROW;
        static constexpr const size_t MaskRadius = 0;
        static constexpr const double MaskErrorRate = 0.0;

        // verbosity & debugging
        static constexpr const bool AnnotateGFF = false;
        static constexpr const bool ReportEffectiveCoverage = false;

        // advanced configuration
        static constexpr const bool UsingFancyChunking = true;
        static constexpr const size_t WindowSpan = 500;
        static constexpr const size_t WindowOverhang = 5;
        static constexpr const bool SkipUnrecognizedContigs = false;
        static constexpr const bool ComputeConfidence = true;
        static constexpr const size_t MaxIterations = 40;
        static constexpr const size_t MutationSeparation = 10;
        static constexpr const size_t MutationNeighborhood = 20;
        static constexpr const float ReadStumpinessThreshold = 0.1;
        static constexpr const SortingStrategy Strategy =
            SortingStrategy::LONGEST_AND_STRAND_BALANCED;
        static constexpr const bool Diploid = false;
        static constexpr const size_t MaxPoaCoverage = 11;
        static constexpr const size_t MinPoaCoverage = 3;
        static constexpr const bool PolishDiploid = true;

        // skipUnrecognizedContigs
    };

    // input files
    std::string inputFilename;
    std::string referenceFilename;

    // output files
    std::string fastaFilename;
    std::string fastqFilename;
    std::string gffFilename;
    std::string vcfFilename;

    // parallelism
    size_t numThreads = Defaults::NumThreads;

    // output settings
    size_t minConfidence = Defaults::MinConfidence;
    size_t minCoverage = Defaults::MinCoverage;
    NoCallStyle noCallStyle = Defaults::NoCall;

    // read selection filters
    size_t maxCoverage = Defaults::MaxCoverage;
    double minAccuracy = Defaults::MinAccuracy;
    float minHqRegionSnr = Defaults::MinHqRegionSnr;
    uint8_t minMapQV = Defaults::MinMapQV;
    float minReadScore = Defaults::MinReadScore;
    double minZScore = Defaults::MinZScore;
    std::vector<std::pair<int16_t, int16_t>> barcodes;  // implement me

    // algorithm and parameters
    double maskErrorRate = Defaults::MaskErrorRate;
    size_t maskRadius = Defaults::MaskRadius;
    ConsensusMode mode = Defaults::Mode;

    // diagnostics
    bool annotateGFF = Defaults::AnnotateGFF;
    std::vector<std::string> dumpEvidence;
    std::string evidenceDirectory;
    bool reportEffectiveCoverage = Defaults::ReportEffectiveCoverage;

    // advanced parameters
    bool computeConfidence = Defaults::ComputeConfidence;
    bool diploid = Defaults::Diploid;
    size_t maxIterations = Defaults::MaxIterations;
    size_t maxPoaCoverage = Defaults::MaxPoaCoverage;
    size_t minPoaCoverage = Defaults::MinPoaCoverage;
    size_t mutationNeighborhood = Defaults::MutationNeighborhood;
    size_t mutationSeparation = Defaults::MutationSeparation;
    bool polishDiploid = Defaults::PolishDiploid;
    float readStumpinessThreshold = Defaults::ReadStumpinessThreshold;
    bool skipUnrecognizedContigs = Defaults::SkipUnrecognizedContigs;  // implement me
    SortingStrategy sortStrategy = Defaults::Strategy;
    bool usingFancyChunking = Defaults::UsingFancyChunking;
    size_t windowSpan = Defaults::WindowSpan;
    size_t windowOverhang = Defaults::WindowOverhang;

    std::vector<ReferenceWindow> filterWindows;
    std::string commandLine;

    // settings -> CLI
    static PacBio::CLI::Interface CreateInterface();

    Settings() = default;

    // CLI -> settings
    Settings(const PacBio::CLI::Results& args);
};

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
