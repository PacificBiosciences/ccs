// Author: Derek Barnett

#pragma once

#include <vector>

#include <pacbio/genomicconsensus/experimental/Consensus.h>
#include <pacbio/genomicconsensus/experimental/Variant.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

///
/// \brief The WindowResult struct
///
struct WindowResult
{
    Consensus css;
    std::vector<Variant> variants;
};

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
