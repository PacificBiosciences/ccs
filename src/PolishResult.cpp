// Author: Lance Hepler

#include <pacbio/consensus/PolishResult.h>

namespace PacBio {
namespace Consensus {

PolishResult operator+(const PolishResult& lhs, const PolishResult& rhs)
{
    PolishResult result;
    result.hasConverged = lhs.hasConverged && rhs.hasConverged;
    result.mutationsTested = lhs.mutationsTested + rhs.mutationsTested;
    result.mutationsApplied = lhs.mutationsApplied + rhs.mutationsApplied;
    result.maxAlphaPopulated.insert(result.maxAlphaPopulated.end(), lhs.maxAlphaPopulated.begin(),
                                    lhs.maxAlphaPopulated.end());
    result.maxBetaPopulated.insert(result.maxBetaPopulated.end(), lhs.maxBetaPopulated.begin(),
                                   lhs.maxBetaPopulated.end());
    result.maxNumFlipFlops.insert(result.maxNumFlipFlops.end(), lhs.maxNumFlipFlops.begin(),
                                  lhs.maxNumFlipFlops.end());
    result.maxAlphaPopulated.insert(result.maxAlphaPopulated.end(), rhs.maxAlphaPopulated.begin(),
                                    rhs.maxAlphaPopulated.end());
    result.maxBetaPopulated.insert(result.maxBetaPopulated.end(), rhs.maxBetaPopulated.begin(),
                                   rhs.maxBetaPopulated.end());
    result.maxNumFlipFlops.insert(result.maxNumFlipFlops.end(), rhs.maxNumFlipFlops.begin(),
                                  rhs.maxNumFlipFlops.end());
    return result;
}
}
}  // ::PacBio::Consensus
