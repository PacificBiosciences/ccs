// Author: Derek Barnett

#pragma once

#include <vector>

#include <pacbio/genomicconsensus/experimental/Settings.h>
#include <pacbio/genomicconsensus/experimental/WorkChunk.h>

namespace PacBio {

namespace CLI {
class Results;
}

namespace GenomicConsensus {
namespace experimental {

struct Workflow
{
    ///
    /// \brief EnumerateChunks
    ///
    /// \param name
    /// \param stride
    /// \param filterWindows
    /// \return
    ///
    static std::vector<WorkChunk> EnumerateChunks(
        const std::string& name, const size_t stride,
        const std::vector<ReferenceWindow>& filterWindows);

    ///
    /// \brief EnumerateWindows
    ///
    /// \param name
    /// \param filterWindows
    /// \return
    ///
    static std::vector<ReferenceWindow> EnumerateWindows(
        const std::string& name, const std::vector<ReferenceWindow>& filterWindows);

    ///
    /// \brief EnumerateWindows
    ///
    /// \param name
    /// \param settings
    /// \return
    ///
    static std::vector<ReferenceWindow> EnumerateWindows(const std::string& name,
                                                         const Settings& settings);

    ///
    /// \brief FancyChunks
    ///
    /// Overloaded for Settings.
    ///
    /// \param name
    /// \param settings
    /// \return
    ///
    static std::vector<WorkChunk> FancyChunks(const std::string& name, const Settings& settings);

    ///
    /// \brief SimpleChunks
    ///
    /// \param name
    /// \param settings
    /// \return
    ///
    static std::vector<WorkChunk> SimpleChunks(const std::string& name, const Settings& settings);

    ///
    /// \brief ReferenceNames
    ///
    /// \param settings
    /// \return
    ///
    static std::vector<std::string> ReferenceNames(const Settings& settings);

    ///
    /// \brief Runner
    ///
    /// Main gcpp application workflow.
    ///
    /// \param results
    /// \return
    ///
    static int Runner(const PacBio::CLI::Results& args);
};

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
