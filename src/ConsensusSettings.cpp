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

#include <pacbio/ccs/ConsensusSettings.h>

namespace PacBio {
namespace CCS {
namespace OptionNames {
using PlainOption = Data::PlainOption;
// clang-format off
const PlainOption MaxLength{
    "max_length",
    { "maxLength" },
    "Maximum Subread Length",
    "Maximum length of subreads to use for generating CCS.",
    CLI::Option::IntType(7000)
};
const PlainOption MinLength{
    "min_length",
    { "minLength" },
    "Minimum Subread Length",
    "Minimum length of subreads to use for generating CCS.",
    CLI::Option::IntType(10)
};
const PlainOption MinPasses{
    "min_passes",
    { "minPasses" },
    "Minimum Number of Passes",
    "Minimum number of subreads required to generate CCS.",
    CLI::Option::IntType(3)
};
const PlainOption MinPredictedAccuracy{
    "min_predicted_accuracy",
    { "minPredictedAccuracy" },
    "Minimum Predicted Accuracy",
    "Minimum predicted accuracy in [0, 1].",
    CLI::Option::FloatType(0.9)
};
const PlainOption MinZScore{
    "min_zscore",
    { "minZScore" },
    "Minimum Z Score",
    "Minimum z-score to use a subread. NaN disables this filter.",
    CLI::Option::FloatType(-3.5)
};
const PlainOption MaxDropFraction{
    "max_drop_fraction",
    { "maxDropFraction" },
    "Maximum Dropped Fraction",
    "Maximum fraction of subreads that can be dropped before giving up.",
    CLI::Option::FloatType(0.34)
};
const PlainOption NoPolish{
    "no_polish",
    { "noPolish" },
    "No Polish CCS",
    "Only output the initial template derived from the POA (faster, less accurate).",
    CLI::Option::BoolType()
};
const PlainOption MinReadScore{
    "min_read_score",
    { "minReadScore" },
    "Minimal Read Score",
    "Minimum read score of input subreads.",
    CLI::Option::FloatType(0.75)
};
const PlainOption MinSnr{
    "min_snr",
    { "minSnr" },
    "Minimum SNR",
    "Minimum SNR of input subreads.",
    CLI::Option::FloatType(3.75)
    // See https://github.com/PacificBiosciences/pbccs/issues/86 for a more
    // detailed discussion of this default.)
};
const PlainOption ByStrand{
    "by_strand",
    { "byStrand" },
    "By Strand CCS",
    "Generate a consensus for each strand.",
    CLI::Option::BoolType()
};
const PlainOption ForceOutput{
    "force",
    { "force" },
    "Force overwrite output",
    "Overwrite OUTPUT file if present.",
    CLI::Option::BoolType()
};
const PlainOption PbIndex{ 
    "pbi", 
    { "pbi" }, 
    "Generate PBI file",
    "Generate a .pbi file for the OUTPUT file.",
    CLI::Option::BoolType()
};
const PlainOption Zmws{
    "zmws",
    { "zmws" },
    "Whitelist ZMWs",
    "Generate CCS for the provided comma-separated holenumber ranges only. Default = all",
    CLI::Option::StringType("")
};
const PlainOption ReportFile{
    "report_file",
    { "reportFile" },
    "Report File Output",
    "Where to write the results report.",
    CLI::Option::StringType("ccs_report.txt")
};
const PlainOption NumThreads{
    "num_threads",
    { "numThreads" },
    "Number of Threads",
    "Number of threads to use, 0 means autodetection.",
    CLI::Option::IntType(0)
};
const PlainOption LogFile{
    "log_file",
    { "logFile" },
    "Log to a File",
    "Log to a file, instead of STDERR.",
    CLI::Option::StringType("")
};
const PlainOption LogLevel{
    "log_level",
    { "logLevel" },
    "Set Log Level",
    "Set log level.",
    CLI::Option::StringType("INFO")
};
const PlainOption RichQVs{
    "rich_qvs",
    { "richQVs" },
    "Emit individual QVs",
    "Emit dq, iq, and sq \"rich\" quality tracks.",
    CLI::Option::BoolType()
};
const PlainOption ModelPath{
    "model_path",
    { "modelPath" },
    "Model(s) Path",
    "Path to a model file or directory containing model files.",
    CLI::Option::StringType("")
};
const PlainOption ModelSpec{
    "model_spec",
    { "modelSpec" },
    "Model Override",
    "Name of chemistry or model to use, overriding default selection.",
    CLI::Option::StringType("")
};
// clang-format on
}  // namespace OptionNames

ConsensusSettings::ConsensusSettings(const PacBio::CLI::Results& options)
    : ByStrand(options[OptionNames::ByStrand])
    , ForceOutput(options[OptionNames::ForceOutput])
    , LogFile(std::forward<std::string>(options[OptionNames::LogFile]))
    , LogLevel(std::forward<std::string>(options[OptionNames::LogLevel]))
    , MaxDropFraction(options[OptionNames::MaxDropFraction])
    , MaxLength(options[OptionNames::MaxLength])
    , MinLength(options[OptionNames::MinLength])
    , MinPasses(options[OptionNames::MinPasses])
    , MinPredictedAccuracy(options[OptionNames::MinPredictedAccuracy])
    , MinReadScore(options[OptionNames::MinReadScore])
    , MinSNR(options[OptionNames::MinSnr])

    , MinZScore(options[OptionNames::MinZScore] == nullptr
                    ? NAN
                    : static_cast<float>(options[OptionNames::MinZScore]))
    , ModelPath(std::forward<std::string>(options[OptionNames::ModelPath]))
    , ModelSpec(std::forward<std::string>(options[OptionNames::ModelSpec]))
    , NoPolish(options[OptionNames::NoPolish])
    , NThreads(ThreadCount(options[OptionNames::NumThreads]))
    , PbIndex(options[OptionNames::PbIndex])
    , ReportFile(std::forward<std::string>(options[OptionNames::ReportFile]))
    , RichQVs(options[OptionNames::RichQVs])
    , WlSpec(std::forward<std::string>(options[OptionNames::Zmws]))
{
}

size_t ConsensusSettings::ThreadCount(int n)
{
    const int m = std::thread::hardware_concurrency();

    if (n < 1) return std::max(1, m + n);

    return std::min(m, n);
}

PacBio::CLI::Interface ConsensusSettings::CreateCLI(const std::string& description,
                                                    const std::string& version)
{
    using Option = PacBio::CLI::Option;
    PacBio::CLI::Interface i{"ccs", description, version};

    i.AlternativeToolContractName("pbccs");

    i.AddHelpOption();     // use built-in help output
    i.AddVersionOption();  // use built-in version output

    // clang-format off
    i.AddPositionalArguments({
        {"input",  "Input file.",  "INPUT"},
        {"output", "Output file.", "OUTPUT"}
    });

    i.AddOptions(
    {
        OptionNames::ForceOutput,
        OptionNames::PbIndex,
        OptionNames::Zmws,
        OptionNames::MaxLength,
        OptionNames::MinLength,
        OptionNames::MinPasses,
        OptionNames::MinPredictedAccuracy,
        OptionNames::MinZScore,
        OptionNames::MaxDropFraction,
        OptionNames::MinSnr,
        OptionNames::MinReadScore,
        OptionNames::ByStrand,
        OptionNames::NoPolish,
        OptionNames::RichQVs,
        OptionNames::ReportFile,
        OptionNames::ModelPath,
        OptionNames::ModelSpec,
        OptionNames::NumThreads,
        OptionNames::LogFile,
        OptionNames::LogLevel
    });

    const std::string id = "pbccs.tasks.ccs";
    CLI::ToolContract::Task tcTask(id);
    tcTask.AddOption(OptionNames::MinSnr);
    tcTask.AddOption(OptionNames::MinReadScore);
    tcTask.AddOption(OptionNames::MaxLength);
    tcTask.AddOption(OptionNames::MinLength);
    tcTask.AddOption(OptionNames::MinPasses);
    tcTask.AddOption(OptionNames::MinPredictedAccuracy);
    tcTask.AddOption(OptionNames::MinZScore);
    tcTask.AddOption(OptionNames::MaxDropFraction);
    tcTask.AddOption(OptionNames::NoPolish);
    tcTask.AddOption(OptionNames::ByStrand);
    tcTask.AddOption(OptionNames::ModelPath);
    tcTask.AddOption(OptionNames::ModelSpec);
    tcTask.AddOption(OptionNames::ReportFile);

    tcTask.InputFileTypes({
        {
            "subread_set",
            "SubreadSet",
            "Subread DataSet or .bam file",
            "PacBio.DataSet.SubreadSet"
        }
    });

    tcTask.OutputFileTypes({
        {
            "bam_output",
            "ConsensusReadSet",
            "Output DataSet XML file",
            "PacBio.DataSet.ConsensusReadSet",
            "ccs"
        }
    });

    CLI::ToolContract::Config tcConfig(tcTask);
    i.EnableToolContract(tcConfig);
    // clang-format on

    return i;
}
}
}  // ::PacBio::CCS