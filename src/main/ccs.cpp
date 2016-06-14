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

#include <boost/algorithm/string/join.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/optional.hpp>

#include <OptionParser.h>

#include <pbbam/BamWriter.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/PbiBuilder.h>
#include <pbbam/PbiFilterQuery.h>
#include <pbbam/ReadGroupInfo.h>

#include <pacbio/ccs/Consensus.h>
#include <pacbio/ccs/ExecUtils.h>
#include <pacbio/ccs/Interval.h>
#include <pacbio/ccs/Logging.h>
#include <pacbio/ccs/ReadId.h>
#include <pacbio/ccs/Utility.h>
#include <pacbio/ccs/Version.h>
#include <pacbio/ccs/Whitelist.h>
#include <pacbio/ccs/WorkQueue.h>

#include <pacbio/consensus/Version.h>

using namespace std;
using namespace PacBio::BAM;
using namespace PacBio::CCS;

using boost::none;
using boost::numeric_cast;
using boost::optional;
using optparse::OptionParser;

// these strings are part of the BAM header, they CANNOT contain newlines
#define DESCRIPTION "Generate circular consensus sequences (ccs) from subreads."

namespace PacBio {
namespace CCS {
namespace OptionNames {
constexpr auto ForceOutput = "force";
constexpr auto PbIndex = "pbi";
constexpr auto Zmws = "zmws";
constexpr auto ReportFile = "reportFile";
constexpr auto NumThreads = "numThreads";
constexpr auto LogFile = "logFile";
constexpr auto LogLevel = "logLevel";
}  // namespace OptionNames
}  // namespace CCS
}  // namespace PacBio

typedef ReadType<ReadId> Subread;
typedef ChunkType<ReadId, Subread> Chunk;
typedef ResultType<ConsensusType> Results;

const auto CircularConsensus = &Consensus<Chunk>;

inline std::string QVsToASCII(const std::vector<int>& qvs)
{
    std::string result;
    result.reserve(qvs.size());

    for (const int qv : qvs) {
        result.push_back(static_cast<char>(std::min(std::max(0, qv), 93) + 33));
    }

    return result;
}

void WriteBamRecords(BamWriter& ccsBam, unique_ptr<PbiBuilder>& ccsPbi, Results& counts,
                     Results&& results)
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

        if (ccs.Strand && *(ccs.Strand) == StrandEnum::FORWARD) name << "/fwd";
        if (ccs.Strand && *(ccs.Strand) == StrandEnum::REVERSE) name << "/rev";

        vector<float> snr = {
            static_cast<float>(ccs.SignalToNoise.A), static_cast<float>(ccs.SignalToNoise.C),
            static_cast<float>(ccs.SignalToNoise.G), static_cast<float>(ccs.SignalToNoise.T)};

        tags["RG"] = MakeReadGroupId(*(ccs.Id.MovieName), "CCS");
        tags["zm"] = static_cast<int32_t>(ccs.Id.HoleNumber);
        tags["np"] = static_cast<int32_t>(ccs.NumPasses);
        tags["rq"] = static_cast<float>(ccs.PredictedAccuracy);
        tags["sn"] = snr;

        // deletion, insertion, and substitution QVs
        tags["dq"] = QVsToASCII(ccs.QVs.DeletionQVs);
        tags["iq"] = QVsToASCII(ccs.QVs.InsertionQVs);
        tags["sq"] = QVsToASCII(ccs.QVs.SubstitutionQVs);

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
        tags["mt"] = static_cast<int32_t>(ccs.MutationsTested);
        tags["ma"] = static_cast<int32_t>(ccs.MutationsApplied);
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
                        unique_ptr<PbiBuilder>&& ccsPbi)
{
    Results counts;
    while (queue.ConsumeWith(WriteBamRecords, ref(*ccsBam), ref(ccsPbi), ref(counts)))
        ;
    return counts;
}

void WriteFastqRecords(ofstream& ccsFastq, Results& counts, Results&& results)
{
    counts += results;
    for (const auto& ccs : results) {
        ccsFastq << '@' << *(ccs.Id.MovieName) << '/' << ccs.Id.HoleNumber << "/ccs";

        if (ccs.Strand && *(ccs.Strand) == StrandEnum::FORWARD) ccsFastq << "/fwd";
        if (ccs.Strand && *(ccs.Strand) == StrandEnum::REVERSE) ccsFastq << "/rev";

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

// TODO(lhepler) move this into ConsensusCore2
bool ValidBaseFeatures(const DataSet& ds)
{
    for (const auto& bam : ds.BamFiles()) {
        for (const auto& rg : bam.Header().ReadGroups()) {
            // P6-C4 and S/P1-C1/beta do not require covariates besides SNR
            if (rg.SequencingChemistry() == "P6-C4" || rg.SequencingChemistry() == "S/P1-C1/beta")
                continue;
            // everything else requires IPD and PulseWidth
            else if (!rg.HasBaseFeature(BaseFeature::IPD) ||
                     !rg.HasBaseFeature(BaseFeature::PULSE_WIDTH))
                return false;
        }
    }
    return true;
}

BamHeader PrepareHeader(const OptionParser& parser, int argc, char** argv, const DataSet& ds)
{
    using boost::algorithm::join;

    ProgramInfo program(parser.prog() + "-" + CCS_VERSION);
    program.Name(parser.prog())
        .CommandLine(parser.prog() + " " + join(vector<string>(argv + 1, argv + argc), " "))
        .Description(DESCRIPTION)
        .Version(CCS_VERSION);

    BamHeader header;
    header.PacBioBamVersion("3.0.1").SortOrder("unknown").Version("1.5").AddProgram(program);

    for (const auto& bam : ds.BamFiles()) {
        for (const auto& rg : bam.Header().ReadGroups()) {
            if (rg.ReadType() != "SUBREAD")
                parser.error("invalid input file, READTYPE must be SUBREAD");

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

size_t ThreadCount(int n)
{
    const int m = thread::hardware_concurrency();

    if (n < 1) return max(1, m + n);

    return min(m, n);
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

int main(int argc, char** argv)
{
    using boost::algorithm::join;
    using boost::make_optional;

    SetColumns();

    // args and options
    //
    //
    // clang-format off
    //   clang messes with the way the arg to version() is formatted..
    auto parser =
        OptionParser()
            .usage("usage: %prog [OPTIONS] INPUT OUTPUT")
            .version("%prog " CCS_VERSION " (commit " CCS_GIT_SHA1 ")"
                     "\nConsensusCore2 " CC2_VERSION " (commit " CC2_GIT_SHA1 ")"
                     "\nCopyright (c) 2014-2015 Pacific Biosciences, Inc.\nLicense: 3-BSD")
            .description(DESCRIPTION
                         "\nAdditional documentation: http://github.com/PacificBiosciences/pbccs");
    // clang-format on
    //
    const vector<string> logLevels = {"TRACE", "DEBUG", "INFO",     "NOTICE",
                                      "WARN",  "ERROR", "CRITICAL", "FATAL"};
    const string em = "--";

    parser.add_option(em + OptionNames::ForceOutput)
        .action("store_true")
        .help("Overwrite OUTPUT file if present.");
    parser.add_option(em + OptionNames::PbIndex)
        .action("store_true")
        .help("Generate a .pbi file for the OUTPUT file.");
    parser.add_option(em + OptionNames::Zmws)
        .help(
            "Generate CCS for the provided comma-separated holenumber ranges only. Default = all");

    ConsensusSettings::AddOptions(&parser);

    parser.add_option(em + OptionNames::ReportFile)
        .set_default("ccs_report.txt")
        .help("Where to write the results report. Default = %default");
    parser.add_option(em + OptionNames::NumThreads)
        .type("int")
        .set_default(0)
        .help("Number of threads to use, 0 means autodetection. Default = %default");
    parser.add_option(em + OptionNames::LogFile).help("Log to a file, instead of STDERR.");
    parser.add_option(em + OptionNames::LogLevel)
        .choices(logLevels.begin(), logLevels.end())
        .set_default("INFO")
        .help("Set log level. Default = %default");

    const auto options = parser.parse_args(argc, argv);
    const auto files = parser.args();

    const ConsensusSettings settings(options);

    const bool forceOutput = options.get(OptionNames::ForceOutput);
    const bool pbIndex = options.get(OptionNames::PbIndex);
    const size_t nThreads = ThreadCount(options.get(OptionNames::NumThreads));
    const size_t chunkSize = 1;

    if (static_cast<int>(options.get(OptionNames::MinPasses)) < 1)
        parser.error("option --minPasses: invalid value: must be >= 1");

    // handle --zmws
    //
    //
    optional<Whitelist> whitelist(none);
    const string wlspec(options.get(OptionNames::Zmws));
    try {
        if (!wlspec.empty()) whitelist = Whitelist(wlspec);
    } catch (...) {
        parser.error("option --zmws: invalid specification: '" + wlspec + "'");
    }

    // input validation
    //
    //
    if (files.size() < 1)
        parser.error("missing INPUT and OUTPUT");
    else if (files.size() < 2)
        parser.error("missing OUTPUT");
    else if (files.size() > 2)
        parser.error("too many arguments for INPUT and OUTPUT");

    // pop first file off the list, is OUTPUT file
    const string inputFile = files.front();
    const string outputFile = files.back();

    // verify input file exists
    if (!FileExists(inputFile)) parser.error("INPUT: file does not exist: '" + inputFile + "'");

    // verify output file does not already exist
    if (FileExists(outputFile) && !forceOutput)
        parser.error("OUTPUT: file already exists: '" + outputFile + "'");

    if (settings.ByStrand && settings.NoPolish)
        parser.error("option --byStrand: incompatible with --noPolish");

    // logging
    //
    //
    ofstream logStream;
    {
        string logLevel(options.get(OptionNames::LogLevel));
        string logFile(options.get(OptionNames::LogFile));

        if (!logFile.empty()) {
            logStream.open(logFile);
            Logging::Logger::Default(new Logging::Logger(logStream, logLevel));
        } else {
            Logging::Logger::Default(new Logging::Logger(cerr, logLevel));
        }
        Logging::InstallSignalHandlers();
    }

    // start processing chunks!
    //
    //
    const auto avail = PacBio::Consensus::SupportedChemistries();

    PBLOG_DEBUG << "Found consensus models for: (" << join(avail, ", ") << ')';

    DataSet ds(inputFile);

    // test that all input chemistries are supported
    {
        set<string> used;
        try {
            used = ds.SequencingChemistries();
        } catch (InvalidSequencingChemistryException& e) {
            PBLOG_FATAL << e.what();
            exit(-1);
        }
        vector<string> unavail;

        set_difference(used.begin(), used.end(), avail.begin(), avail.end(),
                       back_inserter(unavail));

        if (!unavail.empty()) {
            PBLOG_FATAL << "Unsupported chemistries found: \"" << join(unavail, "\", \"") << "\""
                        << ", supported chemistries are: \"" << join(avail, "\", \"") << "\"";
            exit(-1);
        }

        PBLOG_DEBUG << "Using consensus models for: (" << join(used, ", ") << ')';
    }

    if (!ValidBaseFeatures(ds)) {
        PBLOG_FATAL << "Missing base features: IPD or PulseWidth";
        exit(-1);
    }

    const auto filter = PbiFilter::FromDataSet(ds);
    unique_ptr<internal::IQuery> query(nullptr);
    if (filter.IsEmpty())
        query.reset(new EntireFileQuery(ds));
    else
        query.reset(new PbiFilterQuery(filter, ds));

    WorkQueue<Results> workQueue(nThreads);
    future<Results> writer;

    const string outputExt = FileExtension(outputFile);
    if (outputExt == "bam") {
        unique_ptr<BamWriter> ccsBam(
            new BamWriter(outputFile, PrepareHeader(parser, argc, argv, ds)));
        unique_ptr<PbiBuilder> ccsPbi(pbIndex ? new PbiBuilder(outputFile + ".pbi") : nullptr);
        writer = async(launch::async, BamWriterThread, ref(workQueue), move(ccsBam), move(ccsPbi));
    } else if (outputExt == "fastq" || outputExt == "fq")
        writer = async(launch::async, FastqWriterThread, ref(workQueue), ref(outputFile));
    else
        parser.error("OUTPUT: invalid file extension: '" + outputExt + "'");

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
            if (chunk && chunk->size() >= chunkSize) {
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
                chunk->emplace_back(Chunk{ReadId(movieNames[movieName], *holeNumber),
                                          vector<Subread>(), SNR(snr[0], snr[1], snr[2], snr[3]),
                                          read.ReadGroup().SequencingChemistry(), barcodes});
            }
        }

        if (skipZmw) continue;

        // check that barcode matches the previous ones, or else...
        if (barcodes && (!read.HasBarcodes() || !read.HasBarcodeQuality() ||
                         read.BarcodeForward() != get<0>(*barcodes) ||
                         read.BarcodeReverse() != get<1>(*barcodes) ||
                         read.BarcodeQuality() != get<2>(*barcodes))) {
            PBLOG_FATAL << "invalid data: \"bc\" or \"bq\" tag did not agree between subreads!";
            exit(-1);
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

        chunk->back().Reads.emplace_back(
            Subread{ReadId(movieNames[movieName], *holeNumber,
                           Interval(read.QueryStart(), read.QueryEnd())),
                    read.Sequence(), std::move(ipd), std::move(pw), read.LocalContextFlags(),
                    read.ReadAccuracy()});
    }

    // run the remaining tasks
    if (chunk && !chunk->empty()) workQueue.ProduceWith(CircularConsensus, move(chunk), settings);

    // wait for the queue to be done
    workQueue.Finalize();

    // wait for the writer thread and get the results counter
    //   then add in the snr/minPasses counts and write the report
    auto counts = writer.get();
    const string reportFile(options.get(OptionNames::ReportFile));

    if (reportFile == "-")
        WriteResultsReport(cout, counts);
    else {
        ofstream stream(reportFile);
        WriteResultsReport(stream, counts);
    }

    return 0;
}
