// Author: Derek Barnett

#include <pacbio/genomicconsensus/experimental/Workflow.h>

#include <fstream>
#include <future>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <thread>
#include <vector>

#include <pbbam/FastaSequenceQuery.h>
#include <pbbam/IndexedFastaReader.h>
#include <pbcopper/cli/Results.h>
#include <pbcopper/logging/Logging.h>

#include <pacbio/genomicconsensus/experimental/Consensus.h>
#include <pacbio/genomicconsensus/experimental/Filters.h>
#include <pacbio/genomicconsensus/experimental/Intervals.h>
#include <pacbio/genomicconsensus/experimental/Output.h>
#include <pacbio/genomicconsensus/experimental/ReferenceWindow.h>
#include <pacbio/genomicconsensus/experimental/Settings.h>
#include <pacbio/genomicconsensus/experimental/WindowResult.h>
#include <pacbio/genomicconsensus/experimental/WorkChunk.h>
#include <pacbio/parallel/WorkQueue.h>

#include "SettingsOptions.h"

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

namespace {

static void Consumer(PacBio::Parallel::WorkQueue<WindowResult>& queue, const Settings& settings)
{
    auto output = std::make_unique<Output>(settings);

    auto ResultOutput = [&](WindowResult&& result) { output->AddResult(std::move(result)); };
    while (queue.ConsumeWith(ResultOutput))
        ;
}

static WindowResult Producer(const WorkChunk& chunk, const Settings& settings)
{
    PBLOG_INFO << "Processing " << chunk.window;

    // actual work
    // return PacBio::GenomicConsensus::experimental::Process(chunk, settings);

    return WindowResult{Consensus{chunk.window, std::string{}, std::vector<uint8_t>{}},
                        std::vector<Variant>{}};
}

}  // anonymous

std::vector<WorkChunk> Workflow::EnumerateChunks(const std::string& name, const size_t stride,
                                                 const std::vector<ReferenceWindow>& filterWindows)
{
    std::vector<WorkChunk> result;
    const auto windows = EnumerateWindows(name, filterWindows);
    for (const auto& win : windows) {
        const auto intervals = SplitInterval(win.interval, stride);
        for (const auto& interval : intervals)
            result.emplace_back(WorkChunk{{name, interval}, true});
    }
    return result;
}

std::vector<ReferenceWindow> Workflow::EnumerateWindows(
    const std::string& name, const std::vector<ReferenceWindow>& filterWindows)
{
    std::vector<ReferenceWindow> result;
    for (const auto& win : filterWindows) {
        if (win.name == name) result.push_back(win);
    }
    return result;
}

std::vector<ReferenceWindow> Workflow::EnumerateWindows(const std::string& name,
                                                        const Settings& settings)
{
    if (!settings.filterWindows.empty()) return EnumerateWindows(name, settings.filterWindows);

    PacBio::BAM::IndexedFastaReader fasta{settings.referenceFilename};
    const auto length = static_cast<size_t>(fasta.SequenceLength(name));
    return {ReferenceWindow{name, {0, length}}};
}

std::vector<WorkChunk> Workflow::FancyChunks(const std::string& name, const Settings& settings)
{
    std::vector<WorkChunk> result;

    const PacBio::BAM::BamFile bam{settings.inputFilename};
    const PacBio::BAM::PbiRawData index{bam.PacBioIndexFilename()};

    const auto windows = EnumerateWindows(name, settings);
    for (const auto& win : windows) {
        const auto readIntervals = FilteredWindowIntervals(index, win, settings.minMapQV);
        const auto coverageIntervals = CoverageIntervals(win.interval, readIntervals);
        for (const auto& ci : coverageIntervals) {
            const bool hasCoverage = ci.coverage >= settings.minCoverage;
            if (hasCoverage) {
                const ReferenceWindow intervalWin{name, ci.interval};
                const auto chunks = EnumerateChunks(name, settings.windowSpan, {intervalWin});
                for (const auto& chunk : chunks)
                    result.push_back(chunk);
            } else
                result.push_back({ReferenceWindow{name, ci.interval}, false});
        }
    }

    return result;
}

std::vector<std::string> Workflow::ReferenceNames(const Settings& settings)
{
    std::vector<std::string> result;
    std::set<std::string> names;

    if (!settings.filterWindows.empty()) {
        // get names from filter list
        for (const auto& win : settings.filterWindows)
            names.insert(win.name);
    } else {
        // get names from FASTA
        PacBio::BAM::FastaSequenceQuery query{settings.referenceFilename};
        for (const auto& seq : query)
            names.insert(seq.Name());
    }

    for (const auto& name : names)
        result.push_back(name);
    return result;
}

int Workflow::Runner(const CLI::Results& args)
{
    // logging
    //
    // Initialize logging as the very first step. This allows us to redirect
    // incorrect CLI usage to a log file.
    std::ofstream logStream;
    {
        const auto logLevel = args.LogLevel();
        const std::string logFile{args[Options::LogFile].get<decltype(logFile)>()};

        using Logger = PacBio::Logging::Logger;

        Logger* logger;
        if (!logFile.empty()) {
            logStream.open(logFile);
            logger = &Logger::Default(new Logger(logStream, logLevel));
        } else {
            logger = &Logger::Default(new Logger(std::cerr, logLevel));
        }
        PacBio::Logging::InstallSignalHandlers(*logger);
    }

    // setup work queue & output thread
    const Settings settings{args};
    PacBio::Parallel::WorkQueue<WindowResult> workQueue{settings.numThreads};
    std::future<void> writer =
        std::async(std::launch::async, Consumer, std::ref(workQueue), std::ref(settings));

    // main loop: add 'work chunks' to work queue
    const auto referenceNames = ReferenceNames(settings);
    for (const auto& name : referenceNames) {
        const auto chunks = [&name, &settings]() {
            if (settings.usingFancyChunking)
                return FancyChunks(name, settings);
            else
                return SimpleChunks(name, settings);
        }();

        for (const auto& chunk : chunks)
            workQueue.ProduceWith(Producer, chunk, settings);
    }

    // wait for worker/output tasks to finish
    workQueue.Finalize();
    writer.get();

    // any final reporting?

    return EXIT_SUCCESS;
}

std::vector<WorkChunk> Workflow::SimpleChunks(const std::string& name, const Settings& settings)
{
    std::vector<WorkChunk> result;
    const auto windows = EnumerateWindows(name, settings);
    for (const auto& win : windows) {
        const auto intervals = SplitInterval(win.interval, settings.windowSpan);
        for (const auto& interval : intervals)
            result.emplace_back(WorkChunk{{name, interval}, true});
    }
    return result;
}

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
