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

    // window filters
    std::vector<ReferenceWindow> filterWindows;

    // settings -> CLI
    static PacBio::CLI::Interface CreateInterface();

    Settings() = default;

    // CLI -> settings
    Settings(const PacBio::CLI::Results& args);
};

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
