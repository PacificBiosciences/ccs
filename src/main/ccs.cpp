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
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <boost/algorithm/string/join.hpp>
#include <boost/optional.hpp>

#include <OptionParser.h>

#include <pbbam/BamWriter.h>
#include <pbbam/EntireFileQuery.h>

#include <pacbio/ccs/Consensus.h>
#include <pacbio/ccs/ExecUtils.h>
#include <pacbio/ccs/Interval.h>
#include <pacbio/ccs/IntervalTree.h>
#include <pacbio/ccs/Logging.h>
#include <pacbio/ccs/ReadId.h>
#include <pacbio/ccs/WorkQueue.h>

using namespace std;
using namespace PacBio::BAM;
using namespace PacBio::CCS;

using boost::none;
using boost::optional;
using optparse::OptionParser;


#define VERSION "0.0.1"


typedef ReadType<ReadId>           Subread;
typedef ChunkType<ReadId, Subread> Chunk;
typedef ConsensusType<ReadId>      CCS;
typedef ResultType<CCS>            Results;


auto const CircularConsensus = &Consensus<Chunk, CCS>;


void Writer(BamWriter& ccsWriter, Results& counts, Results&& results)
{
    counts += results;

    for (const auto& ccs : results)
    {
        BamRecordImpl record;
        TagCollection tags;
        stringstream  name;

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

        vector<float> snr = {
            static_cast<float>(ccs.SignalToNoise.A),
            static_cast<float>(ccs.SignalToNoise.C),
            static_cast<float>(ccs.SignalToNoise.G),
            static_cast<float>(ccs.SignalToNoise.T) };

        tags["RG"] = MakeReadGroupId(*(ccs.Id.MovieName), "CCS");
        tags["zm"] = static_cast<int32_t>(ccs.Id.HoleNumber);
        tags["np"] = static_cast<int32_t>(ccs.NumPasses);
        tags["rq"] = static_cast<int32_t>(1000 * ccs.PredictedAccuracy);
        tags["sn"] = snr;

        // TODO(lhepler) remove these before release
        tags["ms"] = ccs.ElapsedMilliseconds;
        tags["mt"] = static_cast<int32_t>(ccs.MutationsTested);
        tags["ma"] = static_cast<int32_t>(ccs.MutationsApplied);
        tags["rs"] = ccs.StatusCounts;

        record.Name(name.str())
              .SetSequenceAndQualities(ccs.Sequence, ccs.Qualities)
              .Tags(tags);

        ccsWriter.Write(record);
    }
    ccsWriter.TryFlush();
}


Results WriterThread(WorkQueue<Results>& queue, BamWriter& ccsWriter)
{
    Results counts;
    while (queue.ConsumeWith(Writer, ref(ccsWriter), ref(counts)));
    return counts;
}


BamHeader PrepareHeader(const OptionParser& parser, int argc, char **argv, const vector<string>& files)
{
    using boost::algorithm::join;

    ProgramInfo program(parser.prog() + "-" + VERSION);
    program.Name(parser.prog())
           .CommandLine(parser.prog() + " " + join(vector<string>(argv + 1, argv + argc), " "))
           .Description(parser.description())
           .Version(VERSION);

    BamHeader header;
    header.PacBioBamVersion("3.0b7")
          .SortOrder("unknown")
          .Version("1.5")
          .AddProgram(program);

    for (auto file = ++files.begin(); file != files.end(); ++file)
    {
        BamFile bam(*file);

        for (const auto& rg : bam.Header().ReadGroups())
        {
            if (rg.ReadType() != "SUBREAD")
                parser.error("invalid input file, READTYPE must be SUBREAD");

            ReadGroupInfo readGroup(rg.MovieName(), "CCS");
            readGroup.BindingKit(rg.BindingKit())
                     .SequencingKit(rg.SequencingKit())
                     .BasecallerVersion(rg.BasecallerVersion())
                     .FrameRateHz(rg.FrameRateHz());

            header.AddReadGroup(readGroup);
        }
    }

    return header;
}


size_t ThreadCount(int n)
{
    const int m = thread::hardware_concurrency();

    if (n < 1)
        return max(1, m + n);

    return min(m, n);
}


void WriteResultsReport(ostream& report, const Results& counts)
{
    size_t total = counts.Total();

    report << fixed << setprecision(2);

    report << "Success -- CCS generated," << counts.Success
           << "," << 100.0 * counts.Success / total << '%' << std::endl;

    report << "Failed -- Below SNR threshold," << counts.PoorSNR
           << "," << 100.0 * counts.PoorSNR / total << '%' << std::endl;

    report << "Failed -- No insert regions," << counts.NoSubreads
           << "," << 100.0 * counts.NoSubreads / total << '%' << std::endl;

    report << "Failed -- Insert size too small," << counts.TooShort
           << "," << 100.0 * counts.TooShort / total << '%' << std::endl;

    report << "Failed -- Not enough full passes," << counts.TooFewPasses
           << "," << 100.0 * counts.TooFewPasses / total << '%' << std::endl;

    report << "Failed -- CCS did not converge," << counts.NonConvergent
           << "," << 100.0 * counts.NonConvergent / total << '%' << std::endl;

    report << "Failed -- CCS below minimum predicted accuracy," << counts.PoorQuality
           << "," << 100.0 * counts.PoorQuality / total << '%' << std::endl;
}


int main(int argc, char **argv)
{
    using boost::make_optional;

    SetColumns();

    // args and options
    //
    //
    auto parser = OptionParser()
        .usage("usage: %prog [OPTIONS] OUTPUT FILES...")
        .version(string("%prog ") + VERSION + "\nCopyright (c) 2014 Pacific Biosciences, Inc.\nLicense: 3-BSD")
        .description("A poor man's CCS workflow.");

    vector<string> logLevels = { "TRACE", "DEBUG", "INFO", "NOTICE", "WARN", "ERROR", "CRITICAL", "FATAL" };

    parser.add_option("--report").set_default("-").help("Write results report to a file, instead of STDOUT.");
    parser.add_option("--minSnr").type("float").set_default(4).help("Minimum SNR of input subreads. Default = %default");
    parser.add_option("--minReadScore").type("float").set_default(0.75).help("Minimum read score of input subreads. Default = %default");

    ConsensusSettings::AddOptions(&parser);

    parser.add_option("--zmws").set_default("").help("ZMW holenumbers to generate CCS for, give comma-separated ranges. Default = all");
    parser.add_option("--numThreads").type("int").set_default(0).help("Number of threads to use, 0 means autodetection. Default = %default");
    parser.add_option("--chunkSize").type("int").set_default(5).help("Number of CCS jobs to submit simultaneously. Default = %default");
    parser.add_option("--logFile").help("Log to a file, instead of STDERR.");
    parser.add_option("--logLevel").choices(logLevels.begin(), logLevels.end()).set_default("INFO").help("Set log level. Default = %default");

    const auto options = parser.parse_args(argc, argv);
    const auto files   = parser.args();

    const ConsensusSettings settings(options);

    const float minSnr       = options.get("minSnr");
    const float minReadScore = 1000 * static_cast<float>(options.get("minReadScore"));
    const size_t nThreads    = ThreadCount(options.get("numThreads"));
    const size_t chunkSize   = static_cast<size_t>(options.get("chunkSize"));

    // handle --zmws
    //
    //
    optional<IntervalTree> whitelist(none);
    const string wlspec(options.get("zmws"));
    try
    {
        if (!wlspec.empty())
            whitelist = IntervalTree::FromString(wlspec);
    }
    catch (...)
    {
        parser.error("option --zmws: invalid specification: '" + wlspec + "'");
    }

    // input validation
    //
    //
    if (files.size() < 1)
        parser.error("missing OUTPUT");
    else if (files.size() < 2)
        parser.error("missing FILES...");

    // logging
    //
    //
    fstream logStream;
    {
        string logLevel(options.get("logLevel"));
        string logFile(options.get("logFile"));

        if (!logFile.empty())
        {
            logStream.open(logFile, fstream::out);
            Logging::Logger::Default(new Logging::Logger(logStream, logLevel));
        }
        else
        {
            Logging::Logger::Default(new Logging::Logger(cerr, logLevel));
        }
        Logging::InstallSignalHandlers();
    }

    // start processing chunks!
    //
    //
    BamWriter ccsWriter(files.front(), PrepareHeader(parser, argc, argv, files));
    unique_ptr<vector<Chunk>> chunk(new vector<Chunk>());
    map<string, shared_ptr<string>> movieNames;

    WorkQueue<Results> workQueue(nThreads);
    future<Results> writer = async(launch::async, WriterThread, ref(workQueue), ref(ccsWriter));
    size_t poorSNR = 0, tooFewPasses = 0;

    // skip the first file, it's for output
    for (auto file = ++files.begin(); file != files.end(); ++file)
    {
        BamFile bam(*file);
        EntireFileQuery query(bam);

        // use make_optional here to get around spurious warnings from gcc:
        //   https://gcc.gnu.org/bugzilla/show_bug.cgi?id=47679
        optional<int32_t> holeNumber(none);
        bool skipping = false;

        for (const auto& read : query)
        {
            if (static_cast<float>(read.ReadAccuracy()) < minReadScore)
                continue;

            string movieName = read.MovieName();

            if (movieNames.find(movieName) == movieNames.end())
            {
                movieNames[movieName] = make_shared<string>(movieName);
            }

            if (!holeNumber || *holeNumber != read.HoleNumber())
            {
                if (chunk && !chunk->empty() && chunk->back().Reads.size() < settings.MinPasses)
                {
                    PBLOG_DEBUG << "Skipping ZMW " << chunk->back().Id
                                << ", insufficient number of passes ("
                                << chunk->back().Reads.size() << '<' << settings.MinPasses << ')';
                    tooFewPasses += 1;
                    chunk->pop_back();
                }

                if (chunk && chunk->size() >= chunkSize)
                {
                    // Writer(ccsWriter, Consensus(chunk, minLength, maxPoaCov, minPredAcc));
                    workQueue.ProduceWith(CircularConsensus, move(chunk), settings);
                    chunk.reset(new vector<Chunk>());
                }

                holeNumber = read.HoleNumber();
                auto snr = read.SignalToNoise();
                if (whitelist && !whitelist->Contains(*holeNumber))
                {
                    skipping = true;
                }
                else if (*min_element(snr.begin(), snr.end()) < minSnr)
                {
                    PBLOG_DEBUG << "Skipping ZMW " << movieName << '/' << *holeNumber << ", fails SNR threshold (" << minSnr << ')';
                    poorSNR += 1;
                    skipping = true;
                }
                else
                {
                    skipping = false;
                    chunk->emplace_back(
                        Chunk
                        {
                            ReadId(movieNames[movieName], *holeNumber),
                            vector<Subread>(),
                            SNR(snr)
                        });
                }
            }

            if (!skipping)
            {
                chunk->back().Reads.emplace_back(
                    Subread
                    {
                        ReadId(movieNames[movieName], *holeNumber, Interval(read.QueryStart(), read.QueryEnd())),
                        read.Sequence(),
                        read.InsertionQV(),
                        read.LocalContextFlags()
                    });
            }
        }
    }

    // if the last chunk doesn't have enough passes, skip it
    if (chunk && !chunk->empty() && chunk->back().Reads.size() < settings.MinPasses)
    {
        PBLOG_DEBUG << "Skipping ZMW " << chunk->back().Id
                    << ", insufficient number of passes ("
                    << chunk->back().Reads.size() << '<' << settings.MinPasses << ')';
        tooFewPasses += 1;
        chunk->pop_back();
    }

    // run the remaining tasks
    if (chunk && !chunk->empty())
    {
        // Writer(ccsWriter, Consensus(chunk, minLength, maxPoaCov, minPredAcc));
        workQueue.ProduceWith(CircularConsensus, move(chunk), settings);
    }

    // wait for the queue to be done
    workQueue.Finalize();

    // wait for the writer thread and get the results counter
    //   then add in the snr/minPasses counts and write the report
    auto counts = writer.get();
    counts.PoorSNR += poorSNR;
    counts.TooFewPasses += tooFewPasses;
    const string report(options.get("report"));

    if (report == "-")
    {
        WriteResultsReport(cout, counts);
    }
    else
    {
        fstream stream(report, fstream::out);
        WriteResultsReport(stream, counts);
    }

    return 0;
}
