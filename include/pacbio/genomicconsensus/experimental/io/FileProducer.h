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
