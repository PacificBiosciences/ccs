// Author: David Seifert

#include <fstream>
#include <string>

#include <pbcopper/cli/CLI.h>
#include <pbcopper/logging/Logging.h>

#include <pacbio/UnanimityVersion.h>

#include "VariantCallerSettings.h"

// these strings are part of the BAM header, they CANNOT contain newlines
const std::string DESCRIPTION{
    "Compute genomic consensus and call variants relative to the reference."};
const std::string APPNAME{"variantCaller"};

static int Runner(const PacBio::CLI::Results& args)
{
    // logging
    //
    // Initialize logging as the very first step. This allows us to redirect
    // incorrect CLI usage to a log file.
    std::ofstream logStream;
    {
        const auto logLevel = args.LogLevel();
        const std::string logFile = args["log_file"];

        PacBio::Logging::Logger* logger;
        if (!logFile.empty()) {
            logStream.open(logFile);
            logger =
                &PacBio::Logging::Logger::Default(new PacBio::Logging::Logger(logStream, logLevel));
        } else {
            logger =
                &PacBio::Logging::Logger::Default(new PacBio::Logging::Logger(std::cerr, logLevel));
        }
        InstallSignalHandlers(*logger);
    }

    // Get source args
    const std::vector<std::string> files{args.PositionalArguments()};

    // input validation
    if (files.size() == 0) {
        PBLOG_FATAL << "ERROR: Please provide at least one INPUT and one "
                       "OUTPUT file. See --help for more info about positional "
                       "arguments.";
        std::exit(EXIT_FAILURE);
    }

    const std::string inputFile{files.front()};

    const PacBio::GenomicConsensus::GenomicConsensusSettings settings(args);

    const std::string outputFiles{args["output_filename"].get<decltype(outputFiles)>()};
    const std::string referenceFilename{
        args["reference_filename"].get<decltype(referenceFilename)>()};

    /*
    TODO: Handle dispatching to models here
    */

    return EXIT_SUCCESS;
}

// Entry point
int main(int argc, char* argv[])
{
    const auto version =
        PacBio::UnanimityVersion() + " (commit " + PacBio::UnanimityGitSha1() + ")";
    return PacBio::CLI::Run(
        argc, argv,
        PacBio::GenomicConsensus::GenomicConsensusSettings::CreateCLI(DESCRIPTION, version),
        &Runner);
}
