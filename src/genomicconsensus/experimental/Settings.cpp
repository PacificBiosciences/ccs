// Author: Derek Barnett

#include <pacbio/genomicconsensus/experimental/Settings.h>

#include <cstdlib>

#include <stdexcept>

#include <pacbio/UnanimityVersion.h>
#include <pbcopper/logging/Logging.h>
#include <pbcopper/utility/StringUtils.h>
#include <boost/algorithm/string/predicate.hpp>

#include "SettingsOptions.h"
#include "SettingsToolContract.h"

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {
namespace {  // anonymous

void ParseAlgorithmMode(const PacBio::CLI::Results& args, Settings* const settings)
{
    const std::string algoString{args[Options::Algorithm].get<decltype(algoString)>()};

    if (algoString == "arrow")
        settings->mode = ConsensusMode::ARROW;
    else if (algoString == "plurality")
        settings->mode = ConsensusMode::PLURALITY;
    else if (algoString == "poa")
        settings->mode = ConsensusMode::POA;
    else {
        PBLOG_FATAL << "ERROR: Unrecognized algorithm. See --help for more info.";
        exit(EXIT_FAILURE);
    }
}

void ParseBarcodes(const PacBio::CLI::Results& args, Settings* const settings)
{
    const std::string barcodeArg = args[Options::Barcode].get<decltype(barcodeArg)>();
    if (!barcodeArg.empty()) {
        PBLOG_FATAL << "ERROR: Barcode filtering not yet implemented.";
        exit(EXIT_FAILURE);
    }
}

void ParseDumpEvidence(const PacBio::CLI::Results& args, Settings* const settings)
{
    const std::string evidenceDir = args[Options::EvidenceDirectory].get<decltype(evidenceDir)>();
    const std::string dumpEvidenceTypes =
        args[Options::DumpEvidence].get<decltype(dumpEvidenceTypes)>();
    ;
    if (!dumpEvidenceTypes.empty() || !evidenceDir.empty()) {
        PBLOG_FATAL << "ERROR: Evidence dumping not yet implemented.";
        exit(EXIT_FAILURE);
    }
}

void ParseRequiredFilenames(const PacBio::CLI::Results& args, Settings* const settings)
{
    // BAM file
    const auto& positionalArgs = args.PositionalArguments();
    if (positionalArgs.empty()) {
        PBLOG_FATAL << "ERROR: Input BAM must be provided.";
        exit(EXIT_FAILURE);
    }
    settings->inputFilename = positionalArgs.at(0);

    // reference file
    const std::string refFn = args[Options::ReferenceFilename].get<decltype(refFn)>();
    if (refFn.empty()) {
        PBLOG_FATAL << "ERROR: Input reference must be provided.";
        exit(EXIT_FAILURE);
    }
    settings->referenceFilename = refFn;
}

void ParseFilterWindows(const PacBio::CLI::Results& args, Settings* const settings)
{
    const std::string filterWindowString =
        args[Options::ReferenceWindowsAsString].get<decltype(filterWindowString)>();
    const std::string filterWindowFilename =
        args[Options::ReferenceWindowsFromFile].get<decltype(filterWindowFilename)>();
    if (!filterWindowFilename.empty() || !filterWindowString.empty())
        PBLOG_FATAL << "ERROR: Window filtering not yet implemented.";

    //    std::vector<ReferenceWindow> filterWindows;

    //if options.referenceWindowsAsString is None:
    //    options.referenceWindows = ()
    //elif options.skipUnrecognizedContigs:
    //    # This is a workaround for smrtpipe scatter/gather.
    //    options.referenceWindows = []
    //    for s in options.referenceWindowsAsString.split(","):
    //        try:
    //            win = reference.stringToWindow(s)
    //            options.referenceWindows.append(win)
    //        except Exception:
    //            msg = traceback.format_exc()
    //            logging.debug(msg)
    //            pass
    //else:
    //    options.referenceWindows = map(reference.stringToWindow,
    //                                   options.referenceWindowsAsString.split(","))
    //if options.referenceWindowsFromAlignment:
    //    options.referenceWindows = alnFile.refWindows

    //def stringToWindow(s):
    //    assert isLoaded()
    //    if s is None:
    //        return None
    //    m = re.match("(.*):(.*)-(.*)", s)
    //    if m:
    //        refId    = anyKeyToId(m.group(1))
    //        refStart = int(m.group(2))
    //        refEnd   = min(int(m.group(3)), byName[refId].length)
    //    else:
    //        refId    = anyKeyToId(s)
    //        refStart = 0
    //        refEnd   = byName[refId].length
    //    return (refId, refStart, refEnd)

    //def anyKeyToId(stringKey):
    //    assert isLoaded()
    //    if stringKey in byName:
    //        return byName[stringKey].name
    //    elif stringKey in byPacBioName:
    //        return byPacBioName[stringKey].name
    //    elif stringKey.isdigit():
    //        refId = int(stringKey)
    //        # at this  point, refId can still be the old numeric identifier
    //        return byId[refId].name
    //    else:
    //        raise Exception("Unknown reference name: %s" % stringKey)
}

void ParseNoCallStyle(const PacBio::CLI::Results& args, Settings* const settings)
{
    const std::string noCallString{
        args[Options::NoEvidenceConsensusCall].get<decltype(noCallString)>()};

    if (noCallString == "lowercasereference")
        settings->noCallStyle = NoCallStyle::LOWERCASE_REFERENCE;
    else if (noCallString == "reference")
        settings->noCallStyle = NoCallStyle::REFERENCE;
    else if (noCallString == "nocall")
        settings->noCallStyle = NoCallStyle::NO_CALL;
    else {
        PBLOG_FATAL << "ERROR: 'no evidence consensus call' style: " << noCallString;
        exit(EXIT_FAILURE);
    }
}

void ParseOutputFilenames(const PacBio::CLI::Results& args, Settings* const settings)
{
    const std::string outputFilenameString{
        args[Options::OutputFilenames].get<decltype(outputFilenameString)>()};
    if (outputFilenameString.empty()) PBLOG_WARN << "WARNING: No output files provided.";

    const auto outputFiles = Utility::Split(outputFilenameString, ',');
    for (const auto& fn : outputFiles) {
        if (boost::algorithm::iends_with(fn, ".fasta") || boost::algorithm::iends_with(fn, ".fa")) {
            settings->fastaFilename = fn;
        } else if (boost::algorithm::iends_with(fn, ".fastq") ||
                   boost::algorithm::iends_with(fn, ".fq")) {
            settings->fastqFilename = fn;
        } else if (boost::algorithm::iends_with(fn, ".vcf")) {
            settings->vcfFilename = fn;
        } else if (boost::algorithm::iends_with(fn, ".gff")) {
            settings->gffFilename = fn;
        } else {
            PBLOG_FATAL << "ERROR: Unrecognized extension on output file: " << fn;
            exit(EXIT_FAILURE);
        }
    }
}

void ParseSortStrategy(const PacBio::CLI::Results& args, Settings* const settings)
{
    const std::string strategyString{args[Options::SortStrategy].get<decltype(strategyString)>()};

    if (strategyString == "longest_and_strand_balanced")
        settings->sortStrategy = SortingStrategy::LONGEST_AND_STRAND_BALANCED;
    else if (strategyString == "longest")
        settings->sortStrategy = SortingStrategy::LONGEST;
    else if (strategyString == "spanning")
        settings->sortStrategy = SortingStrategy::SPANNING;
    else if (strategyString == "file_order")
        settings->sortStrategy = SortingStrategy::FILE_ORDER;
    else {
        PBLOG_FATAL << "ERROR: Unrecognized read sorting strategy: " << strategyString;
        exit(EXIT_FAILURE);
    }
}

}  // namespace anonymous

Settings::Settings(const PacBio::CLI::Results& args)
    : numThreads{args[Options::NumThreads]}
    , minConfidence{args[Options::MinConfidence]}
    , minCoverage{args[Options::MinCoverage]}
    , maxCoverage{args[Options::MaxCoverage]}
    , minAccuracy{args[Options::MinAccuracy]}
    , minHqRegionSnr{args[Options::MinSnr]}
    , minMapQV{args[Options::MinMapQV]}
    , minReadScore{args[Options::MinReadScore]}
    , minZScore{args[Options::MinZScore]}
    , maskErrorRate{args[Options::MaskErrorRate]}
    , maskRadius{args[Options::MaskRadius]}
    , annotateGFF{args[Options::AnnotateGFF]}
    , reportEffectiveCoverage{args[Options::ReportEffectiveCoverage]}
    , computeConfidence{!args[Options::FastMode]}
    , diploid{args[Options::Diploid]}
    , maxIterations{args[Options::MaxIterations]}
    , maxPoaCoverage{args[Options::MaxPoaCoverage]}
    , minPoaCoverage{args[Options::MinPoaCoverage]}
    , mutationNeighborhood{args[Options::MutationNeighborhood]}
    , mutationSeparation{args[Options::MutationSeparation]}
    , polishDiploid{args[Options::Diploid]}
    , readStumpinessThreshold{args[Options::ReadStumpinessThreshold]}
    , skipUnrecognizedContigs{args[Options::SkipUnrecognizedContigs]}
    , usingFancyChunking{!args[Options::SimpleChunking]}
    , windowSpan{args[Options::WindowSpan]}
    , windowOverhang{args[Options::WindowOverhang]}
    , commandLine{args.InputCommandLine()}
{
    // "complex" arg parsing
    ParseAlgorithmMode(args, this);
    ParseBarcodes(args, this);
    ParseDumpEvidence(args, this);
    ParseRequiredFilenames(args, this);
    ParseFilterWindows(args, this);
    ParseNoCallStyle(args, this);
    ParseOutputFilenames(args, this);
    ParseSortStrategy(args, this);
}

PacBio::CLI::Interface Settings::CreateInterface()
{
    const std::string appName = "gcpp";
    const std::string description =
        "Compute genomic consensus from alignments and call variants relative to the reference.";
    const auto version =
        PacBio::UnanimityVersion() + " (commit " + PacBio::UnanimityGitSha1() + ")";

    PacBio::CLI::Interface i{appName, description, version};
    i.AddHelpOption();
    i.AddLogLevelOption();
    i.AddVersionOption();
    i.AddGroup("Basic required options", RequiredOptions());
    i.AddGroup("Parallelism", ParallelismOptions());
    i.AddGroup("Output filtering", OutputFilterOptions());
    i.AddGroup("Read selection/filtering", ReadSelectionFilterOptions());
    i.AddGroup("Algorithm and parameter settings", AlgorithmOptions());
    i.AddGroup("Verbosity and debugging", DiagnosticOptions());
    i.AddGroup("Advanced configuration options", AdvancedOptions());
    i.AddPositionalArguments(PositionalArguments());
    i.EnableToolContract(ToolContractConfig());
    return i;
}

// clang-format on

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
