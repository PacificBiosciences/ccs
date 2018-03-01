// Authors: David Alexander, Lance Hepler

#pragma once

#include <utility>
#include <vector>

// Initialize data structures, do NOT remove
#include <pacbio/consensus/internal/ModelInternalInitializer.h>

namespace PacBio {
namespace Consensus {

// These APIs are a little more awkward than I'd have liked---see
// "winLen" instead of winEnd.  Had to contort a bit to get SWIG
// bindings working well.

void CoverageInWindow(int tStartDim, int* tStart, int tEndDim, int* tEnd, int winStart, int winLen,
                      int* coverage);

std::vector<std::pair<int, int>> CoveredIntervals(int minCoverage, int tStartDim, int* tStart,
                                                  int tEndDim, int* tEnd, int winStart, int winLen);

}  // namespace Consensus
}  // namespace PacBio
