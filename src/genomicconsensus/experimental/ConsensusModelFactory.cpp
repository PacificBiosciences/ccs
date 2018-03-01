// Author: Derek Barnett

#include <pacbio/genomicconsensus/experimental/ConsensusModelFactory.h>

#include <memory>
#include <stdexcept>

#include <pacbio/genomicconsensus/experimental/arrow/ArrowModel.h>
#include <pacbio/genomicconsensus/experimental/plurality/PluralityModel.h>
#include <pacbio/genomicconsensus/experimental/poa/PoaModel.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

std::unique_ptr<IConsensusModel> ConsensusModelFactory::Create(const ConsensusMode mode)
{
    switch (mode) {
        case ConsensusMode::ARROW:
            return std::make_unique<ArrowModel>();
        case ConsensusMode::PLURALITY:
            return std::make_unique<PluralityModel>();
        case ConsensusMode::POA:
            return std::make_unique<PoaModel>();
        default:
            throw std::runtime_error("Unknown algorithm mode");
    }
}

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
