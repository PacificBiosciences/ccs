// Author: Derek Barnett

#pragma once

#include <vector>

#include <pacbio/genomicconsensus/experimental/Settings.h>
#include <pacbio/genomicconsensus/experimental/WindowResult.h>
#include <pacbio/genomicconsensus/experimental/WorkChunk.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

///
/// \brief Process
///
/// \param chunk
/// \param settings
///
/// \return
///
WindowResult Process(const WorkChunk& chunk, const Settings& settings);

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
