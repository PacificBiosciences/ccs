// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
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

// Author: Lance Hepler, Armin Töpfer
#pragma once

#include <algorithm>
#include <string>
#include <thread>

#include <pbcopper/cli/CLI.h>

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
    std::string LogLevel;
    double MaxDropFraction;
    size_t MaxLength;
    const size_t MaxPoaCoverage = std::numeric_limits<size_t>::max();
    size_t MinLength;
    size_t MinPasses;
    double MinPredictedAccuracy;
    double MinReadScore;
    double MinSNR;
    double MinZScore;
    std::string ModelPath;
    std::string ModelSpec;
    bool NoPolish;
    size_t NThreads;
    bool PbIndex;
    std::string ReportFile;
    bool RichQVs;
    std::string WlSpec;

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