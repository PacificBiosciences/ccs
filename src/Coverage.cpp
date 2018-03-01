// Author: David Alexander

#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>

#include <pacbio/consensus/Coverage.h>

namespace PacBio {
namespace Consensus {

void CoverageInWindow(int tStartDim, int *tStart, int tEndDim, int *tEnd, int winStart, int winLen,
                      int *coverage)
{
    using std::max;
    using std::min;

    assert(tStartDim == tEndDim);

    int nReads = tStartDim;
    int winEnd = winStart + winLen;
    std::fill_n(coverage, winLen, 0);
    for (int read = 0; read < nReads; read++) {
        int tStart_ = tStart[read];
        int tEnd_ = tEnd[read];
        for (int pos = max(tStart_, winStart); pos < min(tEnd_, winEnd); pos++) {
            coverage[pos - winStart] += 1;
        }
    }
}

static constexpr const int CHUNK_SIZE = 10000;

std::vector<std::pair<int, int>> CoveredIntervals(int minCoverage, int tStartDim, int *tStart,
                                                  int tEndDim, int *tEnd, int winStart, int winLen)
{
    assert(tStartDim == tEndDim);
    // assert(isSorted(tStart));  // find out how to get this ... it's C++11

    // Approach: divide into chunks, find coverage in each chunk,
    // then scan for covered windows ... careful to anneal windows
    // spanning chunk boundaries.  We also rely on the sortedness of the
    // tStart to restrict our attention to

    int winEnd = winStart + winLen;
    int coverage[CHUNK_SIZE];
    int currentIntervalStart = -1;
    std::vector<std::pair<int, int>> intervals;

    int startRowInChunk = 0;
    for (int chunkStart = winStart; chunkStart < winEnd; chunkStart += CHUNK_SIZE) {
        int chunkEnd = std::min(chunkStart + CHUNK_SIZE, winEnd);

        // We compute a conservative guess of the rows that are involved in
        // this chunk.  Not every row in the range [startRowInChunk, endRowInChunk)
        // actually overlaps the chunk, but no rows not in that range intersect the chunk.
        // startRowInChunk is computed by scanning from where it was for the last chunk.
        // This is the best we can do within pulling in the nBackRead stuff.
        int endRowInChunk = std::lower_bound(tStart, tStart + tStartDim, chunkEnd) - tStart;
        for (; ((startRowInChunk < endRowInChunk) && (tEnd[startRowInChunk] < chunkStart));
             startRowInChunk++)
            ;

        CoverageInWindow((endRowInChunk - startRowInChunk), tStart + startRowInChunk,
                         (endRowInChunk - startRowInChunk), tEnd + startRowInChunk, chunkStart,
                         CHUNK_SIZE, coverage);
        int j = 0;
        while (j < (chunkEnd - chunkStart)) {
            if (coverage[j] >= minCoverage) {
                if (currentIntervalStart == -1) {
                    currentIntervalStart = chunkStart + j;
                }
            } else {
                if (currentIntervalStart != -1) {
                    intervals.emplace_back(currentIntervalStart, chunkStart + j);
                    currentIntervalStart = -1;
                }
            }
            j++;
        }
    }
    if (currentIntervalStart != -1) {
        intervals.emplace_back(currentIntervalStart, winEnd);
    }
    return intervals;
}

}  // namespace Consensus
}  // namespace PacBio
