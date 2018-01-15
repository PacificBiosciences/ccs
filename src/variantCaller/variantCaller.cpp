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

// Author: David Seifert

#include <fstream>
#include <string>

#include <pbcopper/cli/CLI.h>
#include <pbcopper/logging/Logging.h>

#include <pacbio/UnanimityVersion.h>

#include "VariantCallerSettings.h"

using namespace std::literals::string_literals;  // for std::operator ""s

// these strings are part of the BAM header, they CANNOT contain newlines
const auto DESCRIPTION{"Compute genomic consensus and call variants relative to the reference."s};
const auto APPNAME{"variantCaller"s};

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
