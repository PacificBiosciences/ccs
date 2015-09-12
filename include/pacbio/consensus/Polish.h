
#pragma once

#include <pacbio/consensus/Mutation.h>

namespace PacBio {
namespace Consensus {

// forward declaration
class AbstractIntegrator;

struct PolishConfig
{
    size_t MaximumIterations;
    size_t MutationSeparation;
    size_t MutationNeighborhood;
};

bool Polish(AbstractIntegrator* ai, const PolishConfig& cfg);

} // namespace Consensus
} // namespace PacBio
