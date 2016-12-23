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

// Author: Armin Töpfer

#include <thread>

#include <pacbio/Version.h>
#include <pacbio/data/PlainOption.h>
#include <boost/algorithm/string.hpp>

#include <pacbio/juliet/JulietSettings.h>
#include <pacbio/juliet/TargetConfig.h>

namespace PacBio {
namespace Juliet {
namespace OptionNames {
using PlainOption = Data::PlainOption;
// clang-format off
const PlainOption Region{
    "region",
    { "region", "r"},
    "Region of Interest",
    "Genomic region of interest, reads will be clipped to that region, empty means all reads.",
    CLI::Option::StringType("2253-5096")
};
const PlainOption Output{
    "output",
    { "output", "o"},
    "Output Prefix",
    "Output prefix for generated files [Default: Input file prefix].",
    CLI::Option::StringType("")
};
const PlainOption PValueThreshold{
    "p_value_threshold",
    { "p-value-threshold", "d" },
    "P-Value Threshold",
    "P-value threshold to call SNV.",
    CLI::Option::FloatType(0.01)
};
const PlainOption DRMOnly{
    "only_known_drms",
    { "drm-only", "k" },
    "Only Known DRMs",
    "Only report known DRM positions.",
    CLI::Option::BoolType()
};
const PlainOption Mode{
    "Execution mode",
    { "mode", "m" },
    "Execution mode",
    "Execution mode: amino, base, phasing, or error",
    CLI::Option::StringType("amino")
};
const PlainOption ErrorModel{
    "Error model",
    { "error", "e" },
    "Error model",
    "Error model: FLEA_RQ95 or FLEA_RQ99",
    CLI::Option::StringType("FLEA_RQ99")
};
const PlainOption TargetConfig{
    "Target config",
    { "config", "c" },
    "Target config",
    "Path to the JSON target config, containing regions of interest, the JSON string itself, or a predefined config tag like <HIV>",
    CLI::Option::StringType("<HIV>")
};
// clang-format on
}  // namespace OptionNames

JulietSettings::JulietSettings(const PacBio::CLI::Results& options)
    : InputFiles(options.PositionalArguments())
    , OutputPrefix(std::forward<std::string>(options[OptionNames::Output]))
    , TargetConfigUser(std::forward<std::string>(options[OptionNames::TargetConfig]))
    , PValueThreshold(options[OptionNames::PValueThreshold])
    , DRMOnly(options[OptionNames::DRMOnly])
    , Mode(AnalysisModeFromString(options[OptionNames::Mode]))
    , SelectedErrorModel(ErrorModelFromString(options[OptionNames::ErrorModel]))
{
    SplitRegion(options[OptionNames::Region], &RegionStart, &RegionEnd);
}

size_t JulietSettings::ThreadCount(int n)
{
    const int m = std::thread::hardware_concurrency();

    if (n < 1) return std::max(1, m + n);

    return std::min(m, n);
}

void JulietSettings::SplitRegion(const std::string& region, int* start, int* end)
{
    if (region.compare("") != 0) {
        std::vector<std::string> splitVec;
        boost::split(splitVec, region, boost::is_any_of("-"));
        *start = stoi(splitVec[0]) - 1;
        *end = stoi(splitVec[1]) - 1;
        if (*start < 0 || *end < 0) throw std::runtime_error("Indexing is 1-based");
    }
}

AnalysisMode JulietSettings::AnalysisModeFromString(const std::string& input)
{
    std::string s = input;
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    if (s.find("amino") != std::string::npos || s.find("acid") != std::string::npos)
        return AnalysisMode::AMINO;
    else if (s.find("base") != std::string::npos || s.find("nuc") != std::string::npos)
        return AnalysisMode::BASE;
    else if (s.find("phas") != std::string::npos || s.find("hap") != std::string::npos)
        return AnalysisMode::PHASING;
    else if (s.find("error") != std::string::npos)
        return AnalysisMode::ERROR;
    else
        throw std::runtime_error("Unknown mode " + s);
}

PacBio::CLI::Interface JulietSettings::CreateCLI()
{
    using Option = PacBio::CLI::Option;
    PacBio::CLI::Interface i{
        "juliet",
        "Juliet, minimal minor variant calling software.\nAttention: Juliet is for research usage "
        "only. Predictions have not been validated.",
        PacBio::UnanimityVersion() + " (commit " + PacBio::UnanimityGitSha1() + ")"};

    i.AddHelpOption();     // use built-in help output
    i.AddVersionOption();  // use built-in version output

    // clang-format off
    i.AddPositionalArguments({
        {"source", "Source BAM or DataSet XML file.", "FILE"}
    });

    i.AddOptions(
    {
        OptionNames::Output,
        OptionNames::Mode,
        OptionNames::ErrorModel,
        OptionNames::Region,
        OptionNames::PValueThreshold,
        OptionNames::DRMOnly,
        OptionNames::TargetConfig
    });
    // clang-format on

    return i;
}
}
}  // ::PacBio::CCS