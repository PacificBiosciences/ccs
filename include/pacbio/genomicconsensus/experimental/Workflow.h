// Copyright (c) 2017, Pacific Biosciences of California, Inc.
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

// Authors: Derek Barnett

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
