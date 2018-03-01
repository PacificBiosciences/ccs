// Author: Derek Barnett

#pragma once

#include <pacbio/genomicconsensus/experimental/WindowResult.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

struct Settings;
struct WorkChunk;

///
/// \brief The IConsensusModel class
///
class IConsensusModel
{
public:
    virtual ~IConsensusModel() = default;

    ///
    /// \brief ProcessChunk
    /// \param chunk
    /// \param settings
    /// \return
    ///
    virtual WindowResult ProcessChunk(const WorkChunk& chunk, const Settings& settings) = 0;

protected:
    IConsensusModel() = default;
};

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
