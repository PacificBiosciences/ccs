// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

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
namespace SparseAlignment {

/// Count the number of seeds in a container.  Additionally,
/// if the MERGESEEDS pre-processing directive is on, adjust
/// the count for the fact that each individual seed may be
/// a composite of multiple smaller seeds.
///
/// \param  seeds  Any container-type of seeds that supports
///                range-based iteration.
template<size_t TSize, typename TContainer>
size_t CountSeeds(const TContainer& seeds)
{
    using namespace seqan;

    size_t count = length(seeds);

#ifdef MERGESEEDS
    for (const auto& seed : seeds)
    {
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
template<size_t TSize>
void FilterSeeds(std::map<size_t, seqan::SeedSet<seqan::Seed<seqan::Simple>>>* seeds,
                 const size_t nBest)
{
    using namespace std;

    // If we already have fewer
    if (seeds->size() <= nBest)
        return;

    // keep a priority queue of the biggest hits,
    // sorted ascendingly. Bump the least value if a new one is bigger.
    priority_queue<size_t, std::vector<size_t>, std::greater<size_t>> best;

    for (const auto& kv : *seeds)
    {
        size_t nSeeds = CountSeeds<TSize>(kv.second);

        if (best.size() < nBest)
        {
            best.push(nSeeds);
        }
        else if (nSeeds > best.top())
        {
            best.pop();
            best.push(nSeeds);
        }
    }

    // Erase all SeedSets with fewer seeds that the smallest
    // item that made it into the queue.
    size_t minSize = best.top();
    for (auto it = seeds->begin(); it != seeds->end(); )
    {
        if (CountSeeds<TSize>(it->second) < minSize)
        {
            it = seeds.erase(it);
        }
        else
        {
            ++it;
        }
    }
}

}  // SparseAlignment
}  // PacBio
