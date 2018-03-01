// Author: Derek Barnett

#pragma once

#include <pacbio/genomicconsensus/experimental/ReferenceWindow.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

///
/// \brief The WorkChunk struct
///
struct WorkChunk
{
    ReferenceWindow window;
    bool hasCoverage;
};

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
