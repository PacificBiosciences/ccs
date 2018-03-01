// Author: David Alexander

#pragma once

// Initialize data structures, do NOT remove
#include <pacbio/consensus/internal/ModelInternalInitializer.h>

namespace PacBio {
namespace Consensus {

enum class MatrixViewConvention
{

    AS_IS,  // View the matrix entries exactly as
            // stored.

    LOGSPACE,  // View matrix entries in logspace;
               // includes column scaling factors.

    LOGPROBABILITY,  // View matrix entries as log-scaled
                     // probabilities.  This entails scaling
                     // per-column and per-row
                     // (counterweights...).

};
}
}
