// Copyright (c) 2017-2018, Pacific Biosciences of California, Inc.
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

// Author: Derek Barnett

#pragma once

#include <cstdio>
#include <string>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

//
// For writers, access sink filenames through this wrapper's tempFilename.
//
// If producer is destructed normally (no live exception), then the temp file
// is renamed to the target name. This provides a clean marker on files that may
// be truncated due to program failure.
//
// For programs that use '-' to indicate writing to stdout, renaming will be
// skipped.
//
struct FileProducer
{
    std::string targetFilename_;
    std::string tempFilename_;

    FileProducer(const std::string& targetFilename)
        : FileProducer{targetFilename, targetFilename + ".tmp"}
    {
    }

    FileProducer(const std::string& targetFilename, const std::string& tempFilename)
        : targetFilename_{targetFilename}, tempFilename_{tempFilename}
    {
    }

    ~FileProducer(void)
    {
        if ((std::current_exception() == nullptr) && (targetFilename_ != "-")) {
            std::rename(tempFilename_.c_str(), targetFilename_.c_str());
        }
    }
};

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
