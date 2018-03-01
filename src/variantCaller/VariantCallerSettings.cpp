// Author: David Seifert

#include "VariantCallerSettings.h"

namespace PacBio {
namespace GenomicConsensus {
namespace OptionNames {

using PlainOption = Data::PlainOption;
// clang-format off

const PlainOption LogFile{
    "log_file",
    { "logFile" },
    "Log to a File",
    "Log to a file, instead of STDERR.",
    CLI::Option::StringType("")
};

// Input/Output files

const PlainOption ReferenceFilename{
    "reference_filename",
    { "r", "referenceFilename" },
    "Reference FASTA",
    "The filename of the reference FASTA file.",
    CLI::Option::StringType("")
};

const PlainOption OutputFilename{
    "output_filename",
    { "o", "outputFilename" },
    "Output files",
    "The output filename(s), as a comma-separated list. Valid output formats "
    "are .fa/.fasta, .fq/.fastq, .gff, .vcf",
    CLI::Option::StringType("")
};

// Parallelism

const PlainOption NumThreads{
    "num_threads",
    { "j", "numWorkers" },
    "Number of Threads",
    "Number of threads to use, 0 means autodetection.",
    CLI::Option::IntType(0)
};

// OutputFiltering

const PlainOption MinConfidence{
    "min_confidence",
    { "q", "minConfidence" },
    "Minimum output confidence",
    "The minimum confidence for a variant call to be output to "
    "variants.{gff,vcf}",
    CLI::Option::IntType(40)
};

const PlainOption MinCoverage{
    "min_coverage",
    { "x", "minCoverage" },
    "Minimum output confidence",
    "The minimum site coverage that must be achieved for variant calls and "
    "consensus to be calculated for a site.",
    CLI::Option::IntType(5)
};

const PlainOption NoEvidenceConsensusCall{
    "no_evidence_consensus_call",
    { "noEvidenceConsensusCall" },
    "Output consensus base",
    "The consensus base that will be output for sites with no effective "
    "coverage. Has to be one of {nocall,reference,lowercasereference}.",
    CLI::Option::StringType("lowercasereference")
};

// ReadSelectionFiltering

const PlainOption Coverage{
    "coverage",
    { "X", "coverage" },
    "Maximum coverage level",
    "A designation of the maximum coverage level to be used for analysis. "
    "Exact interpretation is algorithm-specific.",
    CLI::Option::IntType(100)
};

const PlainOption MinMapQV{
    "min_map_qv",
    { "m", "minMapQV" },
    "Minimum MAPQ",
    "The minimum MapQV for reads that will be used for analysis.",
    CLI::Option::IntType(10)
};

const PlainOption ReferenceWindows{
    "reference_windows",
    { "w", "referenceWindows" },
    "List of reference windows",
    "The window (or multiple comma-delimited windows) of the reference to be "
    "processed, in the format refGroup:refStart-refEnd (default: entire "
    "reference).",
    CLI::Option::StringType("")
};

const PlainOption AlignmentSetRefWindows{
    "alignment_set_ref_windows",
    { "alignmentSetRefWindows" },
    "Load reference windows from file",
    "The window (or multiple comma-delimited windows) of the reference to be "
    "processed, in the format refGroup:refStart-refEnd will be pulled from the "
    "alignment file.",
    CLI::Option::BoolType(false)
};

const PlainOption ReferenceWindowsFile{
    "reference_windows_file",
    { "W", "referenceWindowsFile" },
    "File with list of reference windows",
    "A file containing reference window designations, one per line",
    CLI::Option::StringType("")
};

const PlainOption Barcode{
    "barcode",
    { "barcode" },
    "Barcoded reads to process",
    "Only process reads with the given barcode name.",
    CLI::Option::StringType("")
};

const PlainOption ReadStratum{
    "read_stratum",
    { "readStratum" },
    "Quiver read stratum",
    "A string of the form 'n/N', where n, and N are integers, 0 <= n < N, "
    "designating that the reads are to be deterministically split into N "
    "strata of roughly even size, and stratum n is to be used for variant and "
    "consensus calling. This is mostly useful for Quiver development.",
    CLI::Option::StringType("")
};

const PlainOption MinReadScore{
    "min_read_score",
    { "minReadScore" },
    "Arrow minimum ReadScore",
    "The minimum ReadScore for reads that will be used for analysis "
    "(arrow-only).",
    CLI::Option::FloatType(0.65)
};

const PlainOption MinSnr{
    "min_snr",
    { "minSnr" },
    "Arrow minimum SNR",
    "The minimum acceptable signal-to-noise over all channels for reads that "
    "will be used for analysis (arrow-only).",
    CLI::Option::FloatType(3.75)
};

const PlainOption MinZScore{
    "min_z_score",
    { "minZScore" },
    "Arrow minimum Z-score",
    "The minimum acceptable z-score for reads that will be used for analysis "
    "(arrow-only).",
    CLI::Option::FloatType(-3.5)
};

const PlainOption MinAccuracy{
    "min_accuracy",
    { "minAccuracy" },
    "Arrow minimum accuracy score",
    "The minimum acceptable window-global alignment accuracy for reads that "
    "will be used for the analysis (arrow-only).",
    CLI::Option::FloatType(0.82)
};

// AlgorithmParameterSettings

const PlainOption Algorithm{
    "algorithm",
    { "algorithm" },
    "Used algorithm",
    "The algorithm to use, one of {quiver,arrow,plurality,poa,best}.",
    CLI::Option::StringType("")
};

const PlainOption ParametersFile{
    "parameters_file",
    { "P", "parametersFile" },
    "File with parameter set",
    "Parameter set filename (such as ArrowParameters.json or "
    "QuiverParameters.ini), or directory D such that either D/GenomicConsensus/"
    "QuiverParameters.ini, or D/GenomicConsensus/QuiverParameters.ini, is "
    "found. In the former case, the lexically largest path is chosen.",
    CLI::Option::StringType("")
};

const PlainOption ParametersSpec{
    "parameters_spec",
    { "p", "parametersSpec" },
    "Chemistry model to use",
    "Name of parameter set (chemistry.model) to select from the parameters "
    "file, or just the name of the chemistry, in which case the best available "
    "model is chosen. Default is 'auto', which selects the best parameter set "
    "from the alignment data.",
    CLI::Option::StringType("auto")
};

const PlainOption MaskRadius{
    "mask_radius",
    { "maskRadius" },
    "Mask radius",
    "Radius of window to use when excluding local regions for exceeding "
    "maskMinErrorRate, where 0 disables any filtering (arrow-only).",
    CLI::Option::IntType(3)
};

const PlainOption MaskErrorRate{
    "mask_error_rate",
    { "maskErrorRate" },
    "Maximum allowed error rate before exclusion",
    "Maximum local error rate before the local region defined by maskRadius is "
    "excluded from polishing (arrow-only).",
    CLI::Option::FloatType(0.7)
};

// VerbosityDebuggingProfiling

const PlainOption DumpEvidence{
    "dump_evidence",
    { "d", "dumpEvidence" },
    "Dump variant evidence",
    "Dump evidence relating to variant calling, has to be one of "
    "{variants,all,outliers}.",
    CLI::Option::StringType("")
};

const PlainOption EvidenceDirectory{
    "evidence_directory",
    { "evidenceDirectory" },
    "Directory to dump evidence into",
    "Directory to dump evidence into when enabling --evidenceDirectory.",
    CLI::Option::StringType("")
};

const PlainOption ReportEffectiveCoverage{
    "report_effective_coverage",
    { "reportEffectiveCoverage" },
    "Report effective post-filtering coverage",
    "Additionally record the *post-filtering* coverage at variant sites",
    CLI::Option::BoolType(false)
};

// AdvancedConfiguration

const PlainOption Diploid{
    "diploid",
    { "diploid" },
    "Enable diploid polishing",
    "Enable detection of heterozygous variants (experimental).",
    CLI::Option::BoolType(false)
};

const PlainOption QueueSize{
    "queue_size",
    { "queueSize" },
    "Queue Size",
    "Queue Size",
    CLI::Option::IntType(0)
};

const PlainOption ReferenceChunkSize{
    "reference_chunk_size",
    { "C", "referenceChunkSize" },
    "",
    "",
    CLI::Option::IntType(500)
};

const PlainOption FancyChunking{
    "fancy_chunking",
    { "fancyChunking" },
    "Enable coverage-based chunking",
    "Adaptive reference chunking designed to handle coverage cutouts better",
    CLI::Option::BoolType(true)
};

const PlainOption SimpleChunking{
    "simple_chunking",
    { "simpleChunking" },
    "Disable adaptive reference chunking",
    "Disable adaptive reference chunking",
    CLI::Option::BoolType(true)
};

const PlainOption ReferenceChunkOverlap{
    "reference_chunk_overlap",
    { "referenceChunkOverlap" },
    "Disable adaptive reference chunking",
    "Disable adaptive reference chunking",
    CLI::Option::IntType(10)
};

const PlainOption AutoDisableHdf5ChunkCache{
    "auto_disable_hdf5_chunk_cache",
    { "autoDisableHdf5ChunkCache" },
    "HDF5 chunk cache disabling threshold",
    "Disable the HDF5 chunk cache when the number of datasets in the cmp.h5 "
    "exceeds the given threshold",
    CLI::Option::IntType(500)
};

const PlainOption Aligner{
    "aligner",
    { "a", "aligner" },
    "Quiver variant pairwise alignment algorithm",
    "The pairwise alignment algorithm that will be used to produce variant "
    "calls from the consensus (Quiver only). Has to be one of {affine,simple}.",
    CLI::Option::StringType("affine")
};

const PlainOption RefineDinucleotideRepeats{
    "refine_dinucleotide_repeats",
    { "refineDinucleotideRepeats" },
    "Disable adaptive reference chunking",
    "Require quiver maximum likelihood search to try one less/more repeat copy "
    "in dinucleotide repeats, which seem to be the most frequent cause of "
    "suboptimal convergence (getting trapped in local optimum) (Quiver only)",
    CLI::Option::BoolType(true)
};

const PlainOption NoRefineDinucleotideRepeats{
    "no_refine_dinucleotide_repeats",
    { "noRefineDinucleotideRepeats" },
    "Disable dinucleotide refinement",
    "Disable dinucleotide refinement",
    CLI::Option::BoolType(false)
};

const PlainOption Fast{
    "fast",
    { "fast" },
    "Faster mode",
    "Cut some corners to run faster. Unsupported!",
    CLI::Option::BoolType(false)
};

const PlainOption SkipUnrecognizedContigs{
    "skip_unrecognized_contigs",
    { "skipUnrecognizedContigs" },
    "Ignore reference windows with no coverage",
    "Do not abort when told to process a reference window (via -w/"
    "--referenceWindow[s]) that has no aligned coverage. Outputs emptyish "
    "files if there are no remaining non-degenerate windows. Only intended for "
    "use by smrtpipe scatter/gather.",
    CLI::Option::BoolType(false)
};

// clang-format on

}  // namespace OptionNames

GenomicConsensusSettings::GenomicConsensusSettings(const PacBio::CLI::Results& options)
    : LogFile{options[OptionNames::LogFile].get<decltype(LogFile)>()}
    , LogLevel{options.LogLevel()}
    , ReferenceFilename{options[OptionNames::ReferenceFilename].get<decltype(ReferenceFilename)>()}
    , OutputFilename{options[OptionNames::OutputFilename].get<decltype(OutputFilename)>()}
    , MinConfidence{options[OptionNames::MinConfidence]}
    , MinCoverage{options[OptionNames::MinCoverage]}
    , NoEvidenceConsensusCall{options[OptionNames::NoEvidenceConsensusCall]
                                  .get<decltype(NoEvidenceConsensusCall)>()}
    , Coverage{options[OptionNames::Coverage]}
    , MinMapQV{options[OptionNames::MinMapQV]}
    , ReferenceWindows{options[OptionNames::ReferenceWindows].get<decltype(ReferenceWindows)>()}
    , AlignmentSetRefWindows{options[OptionNames::AlignmentSetRefWindows]}
    , ReferenceWindowsFile{options[OptionNames::ReferenceWindowsFile]
                               .get<decltype(ReferenceWindowsFile)>()}
    , Barcode{options[OptionNames::Barcode].get<decltype(Barcode)>()}
    , ReadStratum{options[OptionNames::ReadStratum].get<decltype(ReadStratum)>()}
    , MinReadScore{options[OptionNames::MinReadScore]}
    , MinSnr{options[OptionNames::MinSnr]}
    , MinZScore{options[OptionNames::MinZScore]}
    , MinAccuracy{options[OptionNames::MinAccuracy]}
    , Algorithm{options[OptionNames::Algorithm].get<decltype(Algorithm)>()}
    , ParametersFile{options[OptionNames::ParametersFile].get<decltype(ParametersFile)>()}
    , ParametersSpec{options[OptionNames::ParametersSpec].get<decltype(ParametersSpec)>()}
    , MaskRadius{options[OptionNames::MaskRadius]}
    , MaskErrorRate{options[OptionNames::MaskErrorRate]}
    , DumpEvidence{options[OptionNames::DumpEvidence].get<decltype(DumpEvidence)>()}
    , EvidenceDirectory{options[OptionNames::EvidenceDirectory].get<decltype(EvidenceDirectory)>()}
    , ReportEffectiveCoverage{options[OptionNames::ReportEffectiveCoverage]}
    , Diploid{options[OptionNames::Diploid]}
    , QueueSize{options[OptionNames::QueueSize]}
    , ReferenceChunkSize{options[OptionNames::ReferenceChunkSize]}
    , FancyChunking{options[OptionNames::FancyChunking]}
    , SimpleChunking{options[OptionNames::SimpleChunking]}
    , ReferenceChunkOverlap{options[OptionNames::ReferenceChunkOverlap]}
    , AutoDisableHdf5ChunkCache{options[OptionNames::AutoDisableHdf5ChunkCache]}
    , Aligner{options[OptionNames::Aligner].get<decltype(Aligner)>()}
    , RefineDinucleotideRepeats{options[OptionNames::RefineDinucleotideRepeats]}
    , NoRefineDinucleotideRepeats{options[OptionNames::NoRefineDinucleotideRepeats]}
    , Fast{options[OptionNames::Fast]}
    , SkipUnrecognizedContigs{options[OptionNames::SkipUnrecognizedContigs]}
{
    // N.B. This is the trick to resolved nthreads from either our
    // option or the "nproc" which has meaning in tool contracts.
    // Derek says he may streamline the API in the future.
    int requestedNThreads;
    if (options.IsFromRTC()) {
        requestedNThreads = options.NumProcessors();
    } else {
        requestedNThreads = options[OptionNames::NumThreads];
    }
    NThreads = ThreadCount(requestedNThreads);
}

size_t GenomicConsensusSettings::ThreadCount(int n)
{
    const int m = std::thread::hardware_concurrency();

    if (n < 1) return std::max(1, m + n);

    return std::min(m, n);
}

PacBio::CLI::Interface GenomicConsensusSettings::CreateCLI(const std::string& description,
                                                           const std::string& version)
{
    using Task = PacBio::CLI::ToolContract::Task;

    PacBio::CLI::Interface i{"variantCaller", description, version};

    i.AlternativeToolContractName("pbgc");

    i.AddHelpOption();      // use built-in help output
    i.AddLogLevelOption();  // use built-in logLevel option
    i.AddVersionOption();   // use built-in version output

    // clang-format off
    i.AddPositionalArguments({
        {"input",  "The input cmp.h5 or BAM alignment file",  "INPUT"}
    });

    i.AddOptions(
    {
        OptionNames::LogFile
    });

    i.AddGroup("Required parameters",
    {
        OptionNames::ReferenceFilename,
        OptionNames::OutputFilename
    });

    i.AddGroup("Parallelism",
    {
        OptionNames::NumThreads
    });

    i.AddGroup("Output filtering",
    {
        OptionNames::MinConfidence,
        OptionNames::MinCoverage,
        OptionNames::NoEvidenceConsensusCall
    });

    i.AddGroup("Read selection/filtering",
    {
        OptionNames::Coverage,
        OptionNames::MinMapQV,
        OptionNames::ReferenceWindows,
        OptionNames::AlignmentSetRefWindows,
        OptionNames::ReferenceWindowsFile,
        OptionNames::Barcode,
        OptionNames::ReadStratum,
        OptionNames::MinReadScore,
        OptionNames::MinSnr,
        OptionNames::MinZScore,
        OptionNames::MinAccuracy
    });

    i.AddGroup("Algorithm and parameter settings",
    {
        OptionNames::Algorithm,
        OptionNames::ParametersFile,
        OptionNames::ParametersSpec,
        OptionNames::MaskRadius,
        OptionNames::MaskErrorRate
    });

    i.AddGroup("Verbosity and debugging/profiling",
    {
        OptionNames::DumpEvidence,
        OptionNames::EvidenceDirectory,
        OptionNames::ReportEffectiveCoverage
    });

    i.AddGroup("Advanced configuration options",
    {
        OptionNames::Diploid,
        OptionNames::QueueSize,
        OptionNames::ReferenceChunkSize,
        OptionNames::FancyChunking,
        OptionNames::SimpleChunking,
        OptionNames::ReferenceChunkOverlap,
        OptionNames::AutoDisableHdf5ChunkCache,
        OptionNames::Aligner,
        OptionNames::RefineDinucleotideRepeats,
        OptionNames::NoRefineDinucleotideRepeats,
        OptionNames::Fast,
        OptionNames::SkipUnrecognizedContigs
    });

    const std::string id = "pbgc.tasks.gc";
    Task tcTask(id);

    tcTask.AddOption(OptionNames::ReferenceFilename);
    tcTask.AddOption(OptionNames::OutputFilename);
    tcTask.AddOption(OptionNames::MinConfidence);
    tcTask.AddOption(OptionNames::MinCoverage);
    tcTask.AddOption(OptionNames::NoEvidenceConsensusCall);
    tcTask.AddOption(OptionNames::Coverage);
    tcTask.AddOption(OptionNames::MinMapQV);
    tcTask.AddOption(OptionNames::ReferenceWindows);
    tcTask.AddOption(OptionNames::AlignmentSetRefWindows);
    tcTask.AddOption(OptionNames::ReferenceWindowsFile);
    tcTask.AddOption(OptionNames::Barcode);
    tcTask.AddOption(OptionNames::ReadStratum);
    tcTask.AddOption(OptionNames::MinReadScore);
    tcTask.AddOption(OptionNames::MinSnr);
    tcTask.AddOption(OptionNames::MinZScore);
    tcTask.AddOption(OptionNames::MinAccuracy);
    tcTask.AddOption(OptionNames::Algorithm);
    tcTask.AddOption(OptionNames::ParametersFile);
    tcTask.AddOption(OptionNames::ParametersSpec);
    tcTask.AddOption(OptionNames::MaskRadius);
    tcTask.AddOption(OptionNames::MaskErrorRate);
    tcTask.AddOption(OptionNames::DumpEvidence);
    tcTask.AddOption(OptionNames::EvidenceDirectory);
    tcTask.AddOption(OptionNames::ReportEffectiveCoverage);
    tcTask.AddOption(OptionNames::Diploid);
    tcTask.AddOption(OptionNames::QueueSize);
    tcTask.AddOption(OptionNames::ReferenceChunkSize);
    tcTask.AddOption(OptionNames::FancyChunking);
    tcTask.AddOption(OptionNames::SimpleChunking);
    tcTask.AddOption(OptionNames::ReferenceChunkOverlap);
    tcTask.AddOption(OptionNames::AutoDisableHdf5ChunkCache);
    tcTask.AddOption(OptionNames::Aligner);
    tcTask.AddOption(OptionNames::RefineDinucleotideRepeats);
    tcTask.AddOption(OptionNames::NoRefineDinucleotideRepeats);
    tcTask.AddOption(OptionNames::Fast);
    tcTask.AddOption(OptionNames::SkipUnrecognizedContigs);
    tcTask.NumProcessors(Task::MAX_NPROC);

    tcTask.InputFileTypes({
        {
            "subread_set",
            "SubreadSet",
            "Aligned Subread DataSet or .bam file",
            "PacBio.DataSet.SubreadSet"
        }
    });

    tcTask.OutputFileTypes({
        {
            "bam_output",
            "Consensus Sequences",
            "Consensus sequences generated by CCS2",
            "PacBio.DataSet.ConsensusReadSet",
            "ccs"
        }
    });

    CLI::ToolContract::Config tcConfig(tcTask);
    i.EnableToolContract(tcConfig);
    // clang-format on

    return i;
}
}  // GC
}  // PacBio
