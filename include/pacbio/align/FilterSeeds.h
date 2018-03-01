// Author: Brett Bowman

#pragma once

#include <functional>
#include <map>
#include <queue>
#include <utility>
#include <vector>

#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/sequence.h>

namespace PacBio {
namespace Align {

/// Count the number of seeds in a container.  Additionally,
/// if the MERGESEEDS pre-processing directive is on, adjust
/// the count for the fact that each individual seed may be
/// a composite of multiple smaller seeds.
///
/// \param  seeds  Any container-type of seeds that supports
///                range-based iteration.
template <size_t TSize, typename TContainer>
size_t CountSeeds(const TContainer& seeds)
{
    using namespace seqan;

    size_t count = length(seeds);

#ifdef MERGESEEDS
    for (const auto& seed : seeds) {
        count += seedSize(seed) - TSize;
    }
#endif

    return count;
}

/// Count the number of seeds in a container.  Additionally,
/// if the MERGESEEDS pre-processing directive is on, adjust
/// the count for the fact that each individual seed may be
/// a composite of multiple smaller seeds.
///
/// \param  seeds  A map of Int<-->SeedSet pairs, linking the index
///                of the reference sequence to the seeds found within.
/// \param  size_t  The number of top hits to retain and report
template <size_t TSize>
void FilterSeeds(std::map<size_t, seqan::SeedSet<seqan::Seed<seqan::Simple>>>* seeds,
                 const size_t nBest)
{
    using namespace std;

    // If we already have fewer
    if (seeds->size() <= nBest) return;

    // keep a priority queue of the biggest hits,
    // sorted ascendingly. Bump the least value if a new one is bigger.
    priority_queue<size_t, std::vector<size_t>, std::greater<size_t>> best;

    for (const auto& kv : *seeds) {
        size_t nSeeds = CountSeeds<TSize>(kv.second);

        if (best.size() < nBest) {
            best.push(nSeeds);
        } else if (nSeeds > best.top()) {
            best.pop();
            best.push(nSeeds);
        }
    }

    // Erase all SeedSets with fewer seeds that the smallest
    // item that made it into the queue.
    size_t minSize = best.top();
    for (auto it = seeds->begin(); it != seeds->end();) {
        if (CountSeeds<TSize>(it->second) < minSize) {
            it = seeds.erase(it);
        } else {
            ++it;
        }
    }
}

}  // Align
}  // PacBio
