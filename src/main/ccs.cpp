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

// Author: Lance Hepler

#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/optional.hpp>

#include <pbcopper/cli/CLI.h>

#include <pbbam/BamWriter.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/PbiBuilder.h>
#include <pbbam/PbiFilterQuery.h>
#include <pbbam/ReadGroupInfo.h>

#include <pbcopper/logging/Logging.h>
#include <pbcopper/utility/FileUtils.h>

#include <pacbio/ccs/Consensus.h>
#include <pacbio/ccs/Whitelist.h>
#include <pacbio/consensus/ModelSelection.h>
#include <pacbio/data/Interval.h>
#include <pacbio/data/ReadId.h>
#include <pacbio/io/Utility.h>
#include <pacbio/parallel/WorkQueue.h>

#include <pacbio/Version.h>

using namespace std;
using namespace PacBio::BAM;
using namespace PacBio::CCS;
using namespace PacBio::Data;
using namespace PacBio::Poa;
using namespace PacBio::Align;
using namespace PacBio::Consensus;
using namespace PacBio::Parallel;
using namespace PacBio::Util;
using namespace PacBio::IO;
using namespace PacBio::Utility;
using namespace PacBio::Logging;

using boost::none;
using boost::numeric_cast;
using boost::optional;

// these strings are part of the BAM header, they CANNOT contain newlines
const string DESCRIPTION = "Generate circular consensus sequences (ccs) from subreads.";
const string APPNAME = "ccs";

typedef ReadType<ReadId> Subread;
typedef ChunkType<ReadId, Subread> Chunk;
typedef ResultType<ConsensusType> Results;

const auto CircularConsensus = &PacBio::CCS::Consensus<Chunk>;

inline string QVsToASCII(const vector<int>& qvs)
{
    string result;
    result.reserve(qvs.size());

    for (const int qv : qvs) {
        result.push_back(static_cast<char>(min(max(0, qv), 93) + 33));
    }

    return result;
}

void WriteBamRecords(BamWriter& ccsBam, unique_ptr<PbiBuilder>& ccsPbi, Results& counts,
                     const bool richQVs, Results&& results)
{
    counts += results;

    for (const auto& ccs : results) {
        BamRecordImpl record;
        TagCollection tags;
        stringstream name;

        // some defaults values
        record.Bin(0)
            .InsertSize(0)
            .MapQuality(255)
            .MatePosition(-1)
            .MateReferenceId(-1)
            .Position(-1)
            .ReferenceId(-1)
            .Flag(0)
            .SetMapped(false);

        name << *(ccs.Id.MovieName) << '/' << ccs.Id.HoleNumber << "/ccs";

        if (ccs.Strand && *(ccs.Strand) == StrandType::FORWARD) name << "/fwd";
        if (ccs.Strand && *(ccs.Strand) == StrandType::REVERSE) name << "/rev";

        vector<float> snr = {
            static_cast<float>(ccs.SignalToNoise.A), static_cast<float>(ccs.SignalToNoise.C),
            static_cast<float>(ccs.SignalToNoise.G), static_cast<float>(ccs.SignalToNoise.T)};

        tags["RG"] = MakeReadGroupId(*(ccs.Id.MovieName), "CCS");
        tags["zm"] = static_cast<int32_t>(ccs.Id.HoleNumber);
        tags["np"] = static_cast<int32_t>(ccs.NumPasses);
        tags["rq"] = static_cast<float>(ccs.PredictedAccuracy);
        tags["sn"] = snr;

        // deletion, insertion, and substitution QVs
        if (richQVs) {
            tags["dq"] = QVsToASCII(ccs.QVs.DeletionQVs);
            tags["iq"] = QVsToASCII(ccs.QVs.InsertionQVs);
            tags["sq"] = QVsToASCII(ccs.QVs.SubstitutionQVs);
        }

        // TODO(lhepler) maybe remove one day
        tags["za"] = static_cast<float>(ccs.AvgZScore);
        vector<float> zScores;
        for (const double z : ccs.ZScores)
            zScores.emplace_back(static_cast<float>(z));
        tags["zs"] = zScores;
        tags["rs"] = ccs.StatusCounts;

        if (ccs.Barcodes) {
            int16_t first, second;
            uint8_t quality;
            tie(first, second, quality) = *ccs.Barcodes;
            vector<uint16_t> bcs{numeric_cast<uint16_t>(first), numeric_cast<uint16_t>(second)};
            tags["bc"] = bcs;
            tags["bq"] = static_cast<int32_t>(quality);
        }

#if DIAGNOSTICS
        tags["ms"] = ccs.ElapsedMilliseconds;
        tags["mt"] = static_cast<int32_t>(ccs.polishResult.mutationsTested);
        tags["ma"] = static_cast<int32_t>(ccs.polishResult.mutationsApplied);
        tags["ap"] = ccs.polishResult.maxAlphaPopulated;
        tags["bp"] = ccs.polishResult.maxBetaPopulated;
        tags["ff"] = ccs.polishResult.maxNumFlipFlops;
#endif

        record.Name(name.str())
            .SetSequenceAndQualities(ccs.Sequence, QVsToASCII(ccs.QVs.Qualities))
            .Tags(tags);

        int64_t offset;
        ccsBam.Write(record, &offset);

        if (ccsPbi) ccsPbi->AddRecord(record, offset);
    }
    ccsBam.TryFlush();
}

Results BamWriterThread(WorkQueue<Results>& queue, unique_ptr<BamWriter>&& ccsBam,
                        unique_ptr<PbiBuilder>&& ccsPbi, const bool richQVs)
{
    Results counts;
    while (queue.ConsumeWith(WriteBamRecords, ref(*ccsBam), ref(ccsPbi), ref(counts), richQVs))
        ;
    return counts;
}

void WriteFastqRecords(ofstream& ccsFastq, Results& counts, Results&& results)
{
    counts += results;
    for (const auto& ccs : results) {
        ccsFastq << '@' << *(ccs.Id.MovieName) << '/' << ccs.Id.HoleNumber << "/ccs";

        if (ccs.Strand && *(ccs.Strand) == StrandType::FORWARD) ccsFastq << "/fwd";
        if (ccs.Strand && *(ccs.Strand) == StrandType::REVERSE) ccsFastq << "/rev";

        ccsFastq << " np:i:" << ccs.NumPasses << " rq:f:" << ccs.PredictedAccuracy;

        if (ccs.Barcodes) {
            int16_t first, second;
            uint8_t quality;
            tie(first, second, quality) = *ccs.Barcodes;
            ccsFastq << " bc:B:S," << first << ',' << second << " bq:i:" << quality;
        }

        ccsFastq << '\n';
        ccsFastq << ccs.Sequence << '\n';
        ccsFastq << "+\n";
        ccsFastq << QVsToASCII(ccs.QVs.Qualities) << '\n';
    }

    ccsFastq.flush();
}

Results FastqWriterThread(WorkQueue<Results>& queue, const string& fname)
{
    ofstream ccsFastq(fname);
    Results counts;
    while (queue.ConsumeWith(WriteFastqRecords, ref(ccsFastq), ref(counts)))
        ;
    return counts;
}

BamHeader PrepareHeader(const string& cmdLine, const DataSet& ds)
{
    using boost::algorithm::join;

    ProgramInfo program(APPNAME + "-" + PacBio::UnanimityVersion());
    program.Name(APPNAME)
        .CommandLine(APPNAME + " " + cmdLine)
        .Description(DESCRIPTION)
        .Version(PacBio::UnanimityVersion());

    BamHeader header;
    header.PacBioBamVersion("3.0.1").SortOrder("unknown").Version("1.5").AddProgram(program);

    for (const auto& bam : ds.BamFiles()) {
        for (const auto& rg : bam.Header().ReadGroups()) {
            if (rg.ReadType() != "SUBREAD") {
                PBLOG_FATAL << "invalid input file, READTYPE must be SUBREAD";
                exit(EXIT_FAILURE);
            }

            ReadGroupInfo readGroup(rg.MovieName(), "CCS");
            readGroup.BindingKit(rg.BindingKit())
                .SequencingKit(rg.SequencingKit())
                .BasecallerVersion(rg.BasecallerVersion())
                .FrameRateHz(rg.FrameRateHz());

            if (rg.HasBarcodeData()) {
                readGroup.BarcodeData(rg.BarcodeFile(), rg.BarcodeHash(), rg.BarcodeCount(),
                                      rg.BarcodeMode(), rg.BarcodeQuality());
            }

            header.AddReadGroup(readGroup);
        }
    }

    return header;
}

void WriteResultsReport(ostream& report, const Results& counts)
{
    size_t total = counts.Total();

    report << fixed << setprecision(2);

    report << "ZMW Yield" << endl;

    report << "Success -- CCS generated," << counts.Success << "," << 100.0 * counts.Success / total
           << '%' << endl;

    report << "Failed -- Below SNR threshold," << counts.PoorSNR << ","
           << 100.0 * counts.PoorSNR / total << '%' << endl;

    report << "Failed -- No usable subreads," << counts.NoSubreads << ","
           << 100.0 * counts.NoSubreads / total << '%' << endl;

    report << "Failed -- Insert size too long," << counts.TooLong << ","
           << 100.0 * counts.TooShort / total << '%' << endl;

    report << "Failed -- Insert size too small," << counts.TooShort << ","
           << 100.0 * counts.TooShort / total << '%' << endl;

    report << "Failed -- Not enough full passes," << counts.TooFewPasses << ","
           << 100.0 * counts.TooFewPasses / total << '%' << endl;

    report << "Failed -- Too many unusable subreads," << counts.TooManyUnusable << ","
           << 100.0 * counts.TooManyUnusable / total << '%' << endl;

    report << "Failed -- CCS did not converge," << counts.NonConvergent << ","
           << 100.0 * counts.NonConvergent / total << '%' << endl;

    report << "Failed -- CCS below minimum predicted accuracy," << counts.PoorQuality << ","
           << 100.0 * counts.PoorQuality / total << '%' << endl;

    report << "Failed -- Unknown error during processing," << counts.ExceptionThrown << ","
           << 100.0 * counts.ExceptionThrown / total << '%' << endl;
    report << endl << endl;

    // Now output the per-subread yield report.
    counts.SubreadCounter.WriteResultsReport(report);
}

static std::vector<ExternalResource> BarcodeSets(const ExternalResources& ext)
{
    std::vector<ExternalResource> output;
    for (const auto& resource : ext) {
        const auto bcs = BarcodeSets(resource.ExternalResources());
        output.insert(output.end(), bcs.begin(), bcs.end());

        if (resource.MetaType() == "PacBio.DataSet.BarcodeSet") output.push_back(resource);
    }
    return output;
}

static int Runner(const PacBio::CLI::Results& args)
{
    // logging
    //
    // Initialize logging as the very first step. This allows us to redirect
    // incorrect CLI usage to a log file.
    ofstream logStream;
    {
        const auto logLevel = args.LogLevel();
        const std::string logFile = args["log_file"];

        Logger* logger;
        if (!logFile.empty()) {
            logStream.open(logFile);
            logger = &Logger::Default(new Logger(logStream, logLevel));
        } else {
            logger = &Logger::Default(new Logger(cerr, logLevel));
        }
        InstallSignalHandlers(*logger);
    }

    using boost::algorithm::join;
    using boost::make_optional;

    // Get source args
    const vector<string> files = args.PositionalArguments();

    // input validation
    if (files.size() != 2) {
        PBLOG_FATAL << "ERROR: Please provide the INPUT and OUTPUT files. See --help for more info "
                       "about positional arguments.";
        exit(EXIT_FAILURE);
    }

    const string inputFile = files.front();
    string outputFile = files.back();

    const ConsensusSettings settings(args);

    // handle --zmws
    //
    //
    optional<Whitelist> whitelist(none);
    const string& wlSpec = settings.WlSpec;
    try {
        if (!wlSpec.empty()) whitelist = Whitelist(wlSpec);
    } catch (const std::exception& e) {
        PBLOG_FATAL << "option --zmws: invalid specification: '" + wlSpec + "'";
        exit(EXIT_FAILURE);
    }

    // verify input file exists
    if (!FileExists(inputFile)) PBLOG_FATAL << "INPUT: file does not exist: '" + inputFile + "'";

    // verify output file does not already exist
    if (FileExists(outputFile) && !settings.ForceOutput) {
        PBLOG_FATAL << "OUTPUT: file already exists: '" + outputFile + "'";
        exit(EXIT_FAILURE);
    }

    if (settings.ByStrand && settings.NoPolish) {
        PBLOG_FATAL << "option --byStrand: incompatible with --noPolish";
        exit(EXIT_FAILURE);
    }

    // load models from file or directory
    //
    //
    {
        const string& modelPath = settings.ModelPath;
        if (!modelPath.empty()) {
            PBLOG_INFO << "Loading model parameters from: '" << modelPath << "'";
            if (!LoadModels(modelPath)) {
                PBLOG_FATAL << "Failed to load models from: " << modelPath;
                exit(EXIT_FAILURE);
            }
        }
    }

    // start processing chunks!
    //
    //
    const auto avail = SupportedChemistries();

    PBLOG_DEBUG << "Found consensus models for: (" << join(avail, ", ") << ')';

    DataSet ds(inputFile);
    const string& modelSpec = settings.ModelSpec;

    // test that all input chemistries are supported
    {
        set<string> used;
        if (!modelSpec.empty()) {
            PBLOG_INFO << "Overriding model selection with: '" << modelSpec << "'";
            if (!(OverrideModel(modelSpec) && used.insert(modelSpec).second)) {
                PBLOG_FATAL << "Failed to find specified model: " << modelSpec;
                exit(EXIT_FAILURE);
            }
        } else {
            try {
                used = ds.SequencingChemistries();
            } catch (InvalidSequencingChemistryException& e) {
                PBLOG_FATAL << e.what();
                exit(EXIT_FAILURE);
            }

            vector<string> unavail;

            set_difference(used.begin(), used.end(), avail.begin(), avail.end(),
                           back_inserter(unavail));

            if (!unavail.empty()) {
                PBLOG_FATAL << "Unsupported chemistries found: (" << join(unavail, ", ") << "), "
                            << "supported chemistries are: (" << join(avail, ", ") << ")";
                exit(EXIT_FAILURE);
            }
        }
        PBLOG_DEBUG << "Using consensus models for: (" << join(used, ", ") << ')';
    }

    if (!ValidBaseFeatures(ds)) {
        PBLOG_FATAL << "Missing base features: IPD or PulseWidth";
        exit(EXIT_FAILURE);
    }

    const auto filter = PbiFilter::FromDataSet(ds);
    unique_ptr<internal::IQuery> query(nullptr);
    if (filter.IsEmpty())
        query.reset(new EntireFileQuery(ds));
    else
        query.reset(new PbiFilterQuery(filter, ds));

    WorkQueue<Results> workQueue(settings.NThreads);
    future<Results> writer;

    // Check if output type is a dataset
    const string outputExt = FileExtension(outputFile);
    bool isXml = outputExt == "xml";
    bool isBam = isXml || outputExt == "bam";

    if (isXml) boost::ireplace_all(outputFile, ".consensusreadset.xml", ".bam");

    if (isBam) {
        unique_ptr<BamWriter> ccsBam(
            new BamWriter(outputFile, PrepareHeader(args.InputCommandLine(), ds)));
        const string pbiFileName = outputFile + ".pbi";
        unique_ptr<PbiBuilder> ccsPbi(new PbiBuilder(pbiFileName));
        writer = async(launch::async, BamWriterThread, ref(workQueue), move(ccsBam), move(ccsPbi),
                       settings.RichQVs);

        // Always generate pbi file
        FileIndex pbi("PacBio.Index.PacBioIndex", pbiFileName);

        if (isXml) {
            // Prepare dataset
            const string metatype = "PacBio.ConsensusReadFile.ConsensusReadBamFile";
            DataSet ccsSet(DataSet::TypeEnum::CONSENSUS_READ);
            ExternalResource resource(metatype, outputFile);

            resource.FileIndices().Add(pbi);

            for (const auto& barcodeSet : BarcodeSets(ds.ExternalResources()))
                resource.ExternalResources().Add(barcodeSet);

            ccsSet.ExternalResources().Add(resource);
            ccsSet.Name(ccsSet.TimeStampedName());

            // File path without .bam suffix
            const auto outputPrefix = outputFile.substr(0, outputFile.size() - 4);
            // Save dataset
            ofstream ccsOut(outputPrefix + ".consensusreadset.xml");
            ccsSet.SaveToStream(ccsOut);
        }
    } else if (outputExt == "fastq" || outputExt == "fq") {
        writer = async(launch::async, FastqWriterThread, ref(workQueue), ref(outputFile));
    } else {
        PBLOG_FATAL << "OUTPUT: invalid file extension: '" + outputExt + "'";
        exit(EXIT_FAILURE);
    }

    unique_ptr<vector<Chunk>> chunk(new vector<Chunk>());
    map<string, shared_ptr<string>> movieNames;
    optional<int32_t> holeNumber(none);
    bool skipZmw = false;
    optional<tuple<int16_t, int16_t, uint8_t>> barcodes(none);

    for (const auto& read : *query) {
        const string movieName = read.MovieName();

        if (movieNames.find(movieName) == movieNames.end())
            movieNames[movieName] = make_shared<string>(movieName);

        // check if we've started a new ZMW
        if (!holeNumber || *holeNumber != read.HoleNumber()) {
            if (chunk && chunk->size() >= settings.ChunkSize) {
                workQueue.ProduceWith(CircularConsensus, move(chunk), settings);
                chunk.reset(new vector<Chunk>());
            }
            holeNumber = read.HoleNumber();
            auto snr = read.SignalToNoise();

            // barcodes
            if (read.HasBarcodes() && read.HasBarcodeQuality()) {
                int16_t first, second;
                uint8_t quality = read.BarcodeQuality();
                tie(first, second) = read.Barcodes();
                barcodes = make_tuple(first, second, quality);
            } else {
                barcodes = none;
            }

            if (whitelist && !whitelist->Contains(movieName, *holeNumber))
                skipZmw = true;
            else {
                skipZmw = false;
                string chem(modelSpec.empty() ? read.ReadGroup().SequencingChemistry() : modelSpec);
                chunk->emplace_back(Chunk{ReadId(movieNames[movieName], *holeNumber),
                                          vector<Subread>(), SNR(snr[0], snr[1], snr[2], snr[3]),
                                          move(chem), barcodes});
            }
        }

        if (skipZmw) continue;

        // check that barcode matches the previous ones, or else...
        if (barcodes && (!read.HasBarcodes() || !read.HasBarcodeQuality() ||
                         read.BarcodeForward() != get<0>(*barcodes) ||
                         read.BarcodeReverse() != get<1>(*barcodes) ||
                         read.BarcodeQuality() != get<2>(*barcodes))) {
            PBLOG_FATAL << "invalid data: \"bc\" or \"bq\" tag did not agree between subreads!";
            exit(EXIT_FAILURE);
        }

        vector<uint8_t> ipd;
        if (read.HasIPD())
            ipd = read.IPD().Encode();
        else
            ipd = vector<uint8_t>(read.Sequence().length(), 0);

        vector<uint8_t> pw;
        if (read.HasPulseWidth())
            pw = read.PulseWidth().Encode();
        else
            pw = vector<uint8_t>(read.Sequence().length(), 0);

        chunk->back().Reads.emplace_back(Subread{
            ReadId(movieNames[movieName], *holeNumber,
                   Interval(read.QueryStart(), read.QueryEnd())),
            read.Sequence(), move(ipd), move(pw), read.LocalContextFlags(), read.ReadAccuracy()});
    }

    // run the remaining tasks
    if (chunk && !chunk->empty()) workQueue.ProduceWith(CircularConsensus, move(chunk), settings);

    // wait for the queue to be done
    workQueue.Finalize();

    // wait for the writer thread and get the results counter
    //   then add in the snr/minPasses counts and write the report
    auto counts = writer.get();
    const string& reportFile = settings.ReportFile;

    if (reportFile == "-")
        WriteResultsReport(cout, counts);
    else {
        ofstream stream(reportFile);
        WriteResultsReport(stream, counts);
    }

    return EXIT_SUCCESS;
}

// Entry point
int main(int argc, char* argv[])
{
    const auto version =
        PacBio::UnanimityVersion() + " (commit " + PacBio::UnanimityGitSha1() + ")";
    return PacBio::CLI::Run(argc, argv, ConsensusSettings::CreateCLI(DESCRIPTION, version),
                            &Runner);
}
