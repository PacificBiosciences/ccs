// Authors: Lance Hepler, Armin Töpfer

#pragma once

#include <algorithm>
#include <string>
#include <thread>

#include <pbcopper/cli/CLI.h>
#include <pbcopper/logging/Logging.h>

#include <pacbio/data/PlainOption.h>

namespace PacBio {
namespace CCS {

/// This class contains all command-line provided arguments and additional
/// constants. Provides a static function to create the CLI pbcopper Interface
/// and the constructor resovlves the CLI::Results automatically.
struct ConsensusSettings
{
    bool ByStrand;
    const size_t ChunkSize = 1;
    bool ForceOutput;
    std::string LogFile;
    Logging::LogLevel LogLevel;
    double MaxDropFraction;
    size_t MaxLength;
    const size_t MaxPoaCoverage = std::numeric_limits<size_t>::max();
    size_t MinLength;
    size_t MinPasses;
    double MinPredictedAccuracy;
    double MinReadScore;
    double MinSNR;
    double MinIdentity;
    double MinZScore;
    std::string ModelPath;
    std::string ModelSpec;
    bool NoPolish;
    size_t PolishRepeats;
    size_t NThreads;
    bool PbIndex;
    std::string ReportFile;
    bool RichQVs;
    std::string WlSpec;
    bool ZmwTimings;

    /// Parses the provided CLI::Results and retrieves a defined set of options.
    ConsensusSettings(const PacBio::CLI::Results& options);

    size_t ThreadCount(int n);

    /// Given the description of the tool and its version, create all
    /// necessary CLI::Options for the ccs executable.
    static PacBio::CLI::Interface CreateCLI(const std::string& description,
                                            const std::string& version);
};
}
}  // ::PacBio::CCS