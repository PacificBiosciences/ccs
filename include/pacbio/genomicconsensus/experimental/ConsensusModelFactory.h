// Author: Derek Barnett

#pragma once

#include <memory>
#include <string>

#include <pacbio/genomicconsensus/experimental/ConsensusMode.h>
#include <pacbio/genomicconsensus/experimental/IConsensusModel.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

///
/// \brief The ConsensusModelFactory struct
///
struct ConsensusModelFactory
{
    ///
    /// \brief Create
    /// \param type
    /// \return
    ///
    static std::unique_ptr<IConsensusModel> Create(const ConsensusMode mode);
};

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
