// Author: David Seifert

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>

#include <pbcopper/cli/CLI.h>
#include <pbcopper/logging/Logging.h>

#include <pacbio/data/PlainOption.h>

namespace PacBio {
namespace GenomicConsensus {

/// This class contains all command-line provided arguments and additional
/// constants. Provides a static function to create the CLI pbcopper Interface
/// and the constructor resovlves the CLI::Results automatically.
struct GenomicConsensusSettings
{
    // Standard
    std::string LogFile;
    Logging::LogLevel LogLevel;

    // Basic required options
    std::string ReferenceFilename;
    std::string OutputFilename;

    // Parallelism
    size_t NThreads;

    // Output filtering
    int32_t MinConfidence;
    int32_t MinCoverage;
    std::string NoEvidenceConsensusCall;

    // Read selection/filtering
    int32_t Coverage;
    int32_t MinMapQV;
    std::string ReferenceWindows;
    bool AlignmentSetRefWindows;
    std::string ReferenceWindowsFile;
    std::string Barcode;
    std::string ReadStratum;
    double MinReadScore;
    double MinSnr;
    double MinZScore;
    double MinAccuracy;

    // Algorithm and parameter settings
    std::string Algorithm;
    std::string ParametersFile;
    std::string ParametersSpec;
    int32_t MaskRadius;
    double MaskErrorRate;

    // Verbosity and debugging/profiling
    std::string DumpEvidence;
    std::string EvidenceDirectory;
    bool ReportEffectiveCoverage;

    // Advanced configuration options
    bool Diploid;
    int32_t QueueSize;
    int32_t ReferenceChunkSize;
    bool FancyChunking;
    bool SimpleChunking;
    int32_t ReferenceChunkOverlap;
    int32_t AutoDisableHdf5ChunkCache;
    std::string Aligner;
    bool RefineDinucleotideRepeats;
    bool NoRefineDinucleotideRepeats;
    bool Fast;
    bool SkipUnrecognizedContigs;

    /// Parses the provided CLI::Results and retrieves a defined set of options.
    GenomicConsensusSettings(const PacBio::CLI::Results& options);

    size_t ThreadCount(int n);

    /// Given the description of the tool and its version, create all
    /// necessary CLI::Options for the ccs executable.
    static PacBio::CLI::Interface CreateCLI(const std::string& description,
                                            const std::string& version);
};
}  // GenomicConsensus
}  // PacBio