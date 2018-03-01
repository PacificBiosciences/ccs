// Author: Derek Barnett

#pragma once

#include <pacbio/data/PlainOption.h>
#include <pacbio/genomicconsensus/experimental/Settings.h>
#include <pbcopper/cli/Option.h>

// clang-format off

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {
namespace Options {

const PacBio::Data::PlainOption ReferenceFilename
{
    "reference_filename",
    {"referenceFilename", "reference", "r"},
    "Reference Filename",
    "The filename of the reference FASTA file.",
    PacBio::CLI::Option::StringType("")
};

const PacBio::Data::PlainOption OutputFilenames
{
    "output_filenames",
    {"outputFilenames", "o"},
    "Output Filenames",
    "The output filename(s), as a comma-separated list. Valid output formats are"
    " .fa/.fasta, .fq/.fastq, .gff, .vcf",
    PacBio::CLI::Option::StringType("")
};

const PacBio::Data::PlainOption NumThreads
{
    "num_threads",
    {"numThreads", "j"},
    "Number of Threads",
    "The number of threads to be used.",
    PacBio::CLI::Option::UIntType(Settings::Defaults::NumThreads)
};

const PacBio::Data::PlainOption MinConfidence
{
    "min_confidence",
    {"minConfidence", "q"},
    "Minimum Confidence",
    "The minimum confidence for a variant call to be output to variants.{gff,vcf}",
    PacBio::CLI::Option::UIntType(Settings::Defaults::MinConfidence)
};

const PacBio::Data::PlainOption MinCoverage
{
    "min_converage",
    {"minCoverage", "x"},
    "Minimum Coverage",
    "The minimum site coverage that must be achieved for variant calls and"
    " consensus to be calculated for a site.",
    PacBio::CLI::Option::UIntType(Settings::Defaults::MinCoverage)
};

const PacBio::Data::PlainOption NoEvidenceConsensusCall
{
    "no_evidence_consensus_call",
    {"noEvidenceConsensusCall"},
    "No Evidence Consensus Call",
    "The consensus base that will be output for sites with no effective coverage.",
    PacBio::CLI::Option::StringType("lowercasereference"),
    {"nocall", "reference", "lowercasereference"}
};

// -------- INPUT FILTER -------- //

const PacBio::Data::PlainOption MaxCoverage
{
    "max_coverage",
    {"coverage", "X"},
    "Maximum Coverage",
    "A designation of the maximum coverage level to be used for analysis. Exact"
    " interpretation is algorithm-specific.",
    PacBio::CLI::Option::UIntType(Settings::Defaults::MaxCoverage)
};

const PacBio::Data::PlainOption MinMapQV
{
    "min_map_qv",
    {"minMapQV", "m"},
    "Minimum MapQV",
    "The minimum MapQV for reads that will be used for analysis.",
    PacBio::CLI::Option::UIntType(Settings::Defaults::MinMapQV)
};

// ** std::vector<ReferenceWindow> filterWindows;
const PacBio::Data::PlainOption ReferenceWindowsAsString
{
    "reference_windows",
    {"referenceWindow", "referenceWindows", "w"},
    "Reference Windows",
    "The window (or multiple comma-delimited windows) of the reference to be"
    " processed, in the format refGroup:refStart-refEnd (default: entire"
    " reference).",
    PacBio::CLI::Option::StringType("")
};

// ** std::vector<ReferenceWindow> filterWindows;
const PacBio::Data::PlainOption ReferenceWindowsFromFile
{
    "reference_windows_from_file",
    {"referenceWindowsFile", "W"},
    "Reference Windows File",
    "A file containing reference window designations, one per line",
    PacBio::CLI::Option::StringType("")
};

const PacBio::Data::PlainOption Barcode
{
    // foo--foo, 0--0

    "barcode",
    {"barcode", "barcodes"},
    "Barcode",
    "Comma-separated list of barcode pairs to analyze, either by name, such as"
    " 'lbc1--lbc1', or by index, such as '0--0'. NOTE: Filtering barcodes by name"
    " requires a barcode file.",
    PacBio::CLI::Option::StringType("")
};

const PacBio::Data::PlainOption BarcodeFile
{
    "barcode_file",
    {"barcodeFile"},
    "Barcode File",
    "Fasta file of the barcode sequences used. NOTE: Only used to find barcode names",
    PacBio::CLI::Option::StringType("")
};

const PacBio::Data::PlainOption MinReadScore
{
    "min_read_score",
    {"minReadScore"},
    "Minimum Read Score",
    "The minimum ReadScore for reads that will be used for analysis (arrow-only).",
    PacBio::CLI::Option::FloatType(Settings::Defaults::MinReadScore)
};

const PacBio::Data::PlainOption MinSnr
{
    "min_hq_region_snr",
    {"minSnr"},
    "Minimum Signal-to-Noise",
    "The minimum acceptable signal-to-noise over all channels for reads that"
    " will be used for analysis (arrow-only).",
    PacBio::CLI::Option::FloatType(Settings::Defaults::MinHqRegionSnr)
};

const PacBio::Data::PlainOption MinZScore
{
    "min_zscore",
    {"minZScore"},
    "Minimum Z-Score",
    "The minimum acceptable z-score for reads that will be used for analysis"
    " (arrow-only).",
    PacBio::CLI::Option::FloatType(Settings::Defaults::MinZScore)
};

const PacBio::Data::PlainOption MinAccuracy
{
    "min_accuracy",
    {"minAccuracy"},
    "Minimum Alignment Accuracy",
    "The minimum acceptable window-global alignment accuracy for reads that will"
    " be used for the analysis (arrow-only).",
    PacBio::CLI::Option::FloatType(Settings::Defaults::MinAccuracy)
};

// -------- ALGORITHM -------- //

const PacBio::Data::PlainOption Algorithm
{
    // GC has a "best" option. Use??

    "algorithm",
    {"algorithm"},
    "Algorithm",
    "The consensus algorithm used.",
    PacBio::CLI::Option::StringType("arrow"),
    {"arrow", "plurality", "poa"}
};

const PacBio::Data::PlainOption ParametersFile
{
    "parameters_file",
    {"parametersFile", "P"},
    "Parameters File",
    "Parameter set filename (such as ArrowParameters.json or QuiverParameters.ini),"
    " or directory D such that either D/*/GenomicConsensus/QuiverParameters.ini,"
    " or D/GenomicConsensus/QuiverParameters.ini, is found.  In the former case,"
    " the lexically largest path is chosen.",
    PacBio::CLI::Option::StringType("")
};

const PacBio::Data::PlainOption ParametersSpec
{
    "parameters_spec",
    {"parametersSpec", "p"},
    "Parameters Spec",
    "Name of parameter set (chemistry.model) to select from the parameters file,"
    " or just the name of the chemistry, in which case the best available model"
    " is chosen.  Default is 'auto', which selects the best parameter set from"
    " the alignment data",
    PacBio::CLI::Option::StringType("auto")
};

const PacBio::Data::PlainOption MaskRadius
{
    "mask_radius",
    {"maskRadius"},
    "Mask Radius",
    "Radius of window to use when excluding local regions for exceeding"
    " maskMinErrorRate, where 0 disables any filtering (arrow-only).",
    PacBio::CLI::Option::UIntType(Settings::Defaults::MaskRadius)
};

const PacBio::Data::PlainOption MaskErrorRate
{
    "mask_error_rate",
    {"maskErrorRate"},
    "Mask Error Rate",
    "Maximum local error rate before the local region defined bymaskRadius is"
    " excluded from polishing (arrow-only).",
    PacBio::CLI::Option::FloatType(Settings::Defaults::MaskErrorRate)
};

const PacBio::Data::PlainOption MaxIterations
{
    "max_iterations",
    {"maxIterations"},
    "Max Iterations",
    "Maximum number of iterations to polish the template.",
    PacBio::CLI::Option::UIntType(Settings::Defaults::MaxIterations)
};

const PacBio::Data::PlainOption MutationSeparation
{
    "mutation_separation",
    {"mutationSeparation"},
    "Mutation Separation",
    "Find the best mutations within a separation window for iterative polishing.",
    PacBio::CLI::Option::UIntType(Settings::Defaults::MutationSeparation)
};

const PacBio::Data::PlainOption MutationNeighborhood
{
    "mutation_neighborhood",
    {"mutationNeighborhood"},
    "Mutation Neighborhood",
    "Find nearby mutations within neighborhood for iterative polishing.",
    PacBio::CLI::Option::UIntType(Settings::Defaults::MutationSeparation)
};

const PacBio::Data::PlainOption ReadStumpinessThreshold
{
    "read_stumpiness_threshold",
    {"readStumpinessThreshold"},
    "Read Stumpinesss Threshold",
    "Filter out reads whose aligned length along a subread is lower than a percentage of its corresponding reference length.",
    PacBio::CLI::Option::FloatType(Settings::Defaults::ReadStumpinessThreshold)
};

const PacBio::Data::PlainOption MaxPoaCoverage
{
    "max_poa_coverage",
    {"maxPoaCoverage"},
    "Maximum POA Coverage",
    "Maximum number of sequences to use for consensus calling.",
    PacBio::CLI::Option::UIntType(Settings::Defaults::MaxPoaCoverage)
};

// -------- DIAGNOSTICS -------- //

const PacBio::Data::PlainOption LogFile
{
    "log_file",
    {"logFile"},
    "Log File",
    "Log to a file, instead of STDERR.",
    CLI::Option::StringType("")
};

const PacBio::Data::PlainOption DumpEvidence
{
    // TODO: how to implement a la argparse(default=None, const="variants", nargs="?", choices=["variants", "all", "outliers"]) ??

    "dump_evidence",
    {"dumpEvidence", "d"},
    "Dump Evidence",
    "Dump evidence data",
    PacBio::CLI::Option::StringType(""),
    {"variants", "all", "outliers"}
};

const PacBio::Data::PlainOption EvidenceDirectory
{
    "evidence_directory",
    {"evidenceDirectory"},
    "Evidence Directory",
    "Directory to dump evidence into.",
    PacBio::CLI::Option::StringType("")
};

const PacBio::Data::PlainOption AnnotateGFF
{
    "annotate_gff",
    {"annotateGFF"},
    "Annotate GFF",
    "Augment GFF variant records with additional information",
    PacBio::CLI::Option::BoolType(Settings::Defaults::AnnotateGFF)
};

const PacBio::Data::PlainOption ReportEffectiveCoverage
{
    "report_effective_coverage",
    {"reportEffectiveCoverage"},
    "Report Effective Coverage",
    "Additionally record the *post-filtering* coverage at variant sites",
    PacBio::CLI::Option::BoolType(Settings::Defaults::ReportEffectiveCoverage)
};

// -------- ADVANCED CONFIG -------- //

const PacBio::Data::PlainOption SortStrategy
{
    //SortingStrategy sortStrategy = Defaults::Strategy;

    "sort_strategy",
    {"sortStrategy"},
    "Read Sorting Strategy",
    "Read sortiing strategy",
    PacBio::CLI::Option::StringType("longest_and_strand_balanced"),
    {"longest_and_strand_balanced", "longest", "spanning", "file_order"}
};

const PacBio::Data::PlainOption Diploid
{
    //
    // TODO: what's the difference between these?
    //
    // bool polishDiploid = Defaults::PolishDiploid;   // ???
    // bool diploid = Defaults::Diploid;
    //

    "diploid",
    {"diploid"},
    "Detect Heterozygous Variants",
    "Enable detection of heterozygous variants (experimental)",
    PacBio::CLI::Option::BoolType(Settings::Defaults::Diploid)
};

const PacBio::Data::PlainOption WindowSpan
{
    // legacy: referenceChunkSize

    "window_span",
    {"referenceChunkSize", "C"},
    "Reference Chunk Size",
    "Size of reference chunks.",
    PacBio::CLI::Option::UIntType(Settings::Defaults::WindowSpan)
};

const PacBio::Data::PlainOption SimpleChunking
{
    "simpleChunking",
    {"simpleChunking"},
    "Simple Chunking",
    "Disable adaptive reference chunking.",
    PacBio::CLI::Option::BoolType(!Settings::Defaults::UsingFancyChunking)
};

const PacBio::Data::PlainOption WindowOverhang
{
    // legacy: referenceChunkOverlap

    "referenceChunkOverlap",
    {"referenceChunkOverlap"},
    "Reference Chunk Overlap",
    "Size of reference chunk overlaps.",
    PacBio::CLI::Option::UIntType(Settings::Defaults::WindowOverhang)
};

const PacBio::Data::PlainOption FastMode
{
    "fast_mode",
    {"fast"},
    "Fast Mode",
    "Cut some corners to run faster.  Unsupported!",
    PacBio::CLI::Option::BoolType(!Settings::Defaults::ComputeConfidence)
};

const PacBio::Data::PlainOption SkipUnrecognizedContigs
{
    "skip_unrecognized_contigs",
    {"skipUnrecognizedContigs"},
    "Skip Unrecognized Contigs",
    "Do not abort when told to process a reference window (via"
    " -w/--referenceWindow[s]) that has no aligned coverage.  Outputs emptyish"
    " files if there are no remaining non-degenerate windows.  Only intended"
    " for use by smrtpipe scatter/gather.",
    PacBio::CLI::Option::BoolType(false)
};

const PacBio::Data::PlainOption MinPoaCoverage
{
    // TODO: difference between this and MinCoverage ??

    "min_poa_coverage",
    {"minPoaCoverage"},
    "Minimum Poa Coverage",
    "Minimum number of reads required within a window to call consensus and"
    " variants using arrow or poa.",
    PacBio::CLI::Option::UIntType(Settings::Defaults::MinPoaCoverage)
};

const PacBio::CLI::PositionalArg InputFilename
{
    "INPUT",
    "The input BAM alignment file",
    "INPUT"
};

}  // namespace Options

inline std::vector<PacBio::CLI::Option> RequiredOptions()
{
    return
    {
        Options::ReferenceFilename,
        Options::OutputFilenames
    };
}

inline std::vector<PacBio::CLI::Option> ParallelismOptions()
{
    return
    {
        Options::NumThreads
    };
}

inline std::vector<PacBio::CLI::Option> OutputFilterOptions()
{
    return
    {
        Options::MinConfidence,
        Options::MinCoverage,
        Options::NoEvidenceConsensusCall
    };
}

inline std::vector<PacBio::CLI::Option> ReadSelectionFilterOptions()
{
    return
    {
        Options::MaxCoverage,
        Options::MinAccuracy,
        Options::MinMapQV,
        Options::MinReadScore,
        Options::MinSnr,
        Options::MinZScore,
        Options::Barcode,
        Options::BarcodeFile,
        Options::ReferenceWindowsAsString,
        Options::ReferenceWindowsFromFile
    };
}

inline std::vector<PacBio::CLI::Option> AlgorithmOptions()
{
    return
    {
        Options::Algorithm,
        Options::MaskRadius,
        Options::MaskErrorRate,
        Options::ParametersFile,
        Options::ParametersSpec,
        Options::MaxIterations,
        Options::MaxPoaCoverage,
        Options::MutationSeparation,
        Options::MutationNeighborhood,
        Options::ReadStumpinessThreshold
    };
}

inline std::vector<PacBio::CLI::Option> DiagnosticOptions()
{
    return
    {
        Options::LogFile,
        Options::DumpEvidence,
        Options::EvidenceDirectory,
        Options::AnnotateGFF,
        Options::ReportEffectiveCoverage
    };
}

inline std::vector<PacBio::CLI::Option> AdvancedOptions()
{
    return
    {
        Options::WindowSpan,
        Options::WindowOverhang,
        Options::SimpleChunking,
        Options::Diploid,
        Options::FastMode,
        Options::SkipUnrecognizedContigs,
        Options::SortStrategy,
        Options::MinPoaCoverage
    };
}

inline std::vector<PacBio::CLI::PositionalArg> PositionalArguments()
{
    return
    {
        Options::InputFilename
    };
}

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio

// clang-format on
