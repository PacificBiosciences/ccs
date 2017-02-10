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

// Author: Lance Hepler

#include <algorithm>
#include <set>

#include <boost/optional.hpp>

#include <pacbio/align/ChainSeeds.h>
#include <pacbio/ccs/ChainSeeds.h>

namespace PacBio {
namespace CCS {
namespace {

// H is query, V is reference

using PacBio::Align::ComputeVisibilityLeft;
using PacBio::Align::Diagonal;
using PacBio::Align::DiagonalCompare;
using PacBio::Align::HVCompare;
using PacBio::Align::IndexCompare;
using PacBio::Align::VHCompare;

using SDPHit = PacBio::Align::SDPHit;
using SDPColumn = PacBio::Align::SDPColumn;

long LinkScore(const PacBio::Align::Seed& a,
               const PacBio::Align::Seed& b,
               const int matchReward)
{
    using namespace std;

    // TODO (dbarnett) : sync this up with PacBio::Align::LinkScore

    long aH = static_cast<long>(a.BeginPositionH());
    long aV = static_cast<long>(a.BeginPositionV());
    long bH = static_cast<long>(b.BeginPositionH());
    long bV = static_cast<long>(b.BeginPositionV());
    long k = min(static_cast<long>(a.Size()), static_cast<long>(b.Size()));

    long diagA = Diagonal(a);
    long diagB = Diagonal(b);
    long fwd = min(aH - bH, aV - bV);
    long indels = abs(diagA - diagB);
    long matches = k - max(0L, k - fwd);
    long mismatches = fwd - matches;

    return matchReward * matches - indels - mismatches;
}

}  // anonymous namespace

std::vector<PacBio::Align::Seed> ChainSeeds(const PacBio::Align::Seeds& seedSet,
                                            const int matchReward)
{
    using namespace std;

    // TODO (dbarnett) : sync this up with PacBio::Align::ChainSeeds

    vector<SDPHit> seeds;
    vector<long> scores(seedSet.size(), 0L);

    // initialize our "SDPHits" vector (has fixed index)
    {
        size_t i = 0;

        for (auto it = seedSet.begin(); it != seedSet.end(); ++it, ++i) {
            seeds.push_back(SDPHit(*it, i));
            scores[i] = (*it).Size();
        }
    }

    // compute visibility left, requires H-sorted seeds
    std::set<SDPHit> sweepSet;
    sort(seeds.begin(), seeds.end(), HVCompare);
    auto visible = ComputeVisibilityLeft(seeds, sweepSet);

    // compute the visibility above, requires V-sorted seeds
    sort(seeds.begin(), seeds.end(), VHCompare);
    auto toRemove = seeds.begin();
    std::set<SDPColumn> colSet;

    long bestChainScore = -numeric_limits<long>::max();
    boost::optional<size_t> bestChainEnd = boost::none;
    vector<boost::optional<size_t>> chainPred(seeds.size(), boost::none);

    auto zScore = [&scores](const SDPHit& seed) {
        return scores[seed.Index] + seed.BeginPositionH() + seed.BeginPositionV();
    };

    for (auto it = seeds.begin(); it != seeds.end();) {
        const size_t row = (*it).BeginPositionV();
        const auto start = it;

        for (; it != seeds.end() && row == (*it).BeginPositionV(); ++it) {
            long bestScore = -numeric_limits<long>::max();
            boost::optional<SDPHit> bestSeed = boost::none;

            // find the previous column and best fragment from it
            {
                SDPColumn col((*it).BeginPositionH());
                auto pred = colSet.lower_bound(col);

                if (pred != colSet.begin()) {
                    pred = prev(pred);
                    long s = scores[pred->Seed->Index] + LinkScore(*it, *(pred->Seed), matchReward);

                    if (s > bestScore) {
                        bestScore = s;
                        bestSeed = pred->Seed;
                    }
                }
            }

            // search visible fragments (above)
            {
                auto visa = sweepSet.lower_bound(*it);

                if (visa != sweepSet.begin()) {
                    visa = prev(visa);
                    long s = scores[visa->Index] + LinkScore(*it, *visa, matchReward);

                    if (s > bestScore) {
                        bestScore = s;
                        bestSeed = *visa;
                    }
                }
            }

            // search visible fragments (left)
            if (auto visl = visible[it->Index]) {
                long s = scores[visl->Index] + LinkScore(*it, *visl, matchReward);

                if (s > bestScore) {
                    bestScore = s;
                    bestSeed = visl;
                }
            }

            if (bestSeed && bestScore > 0) {
                scores[it->Index] = bestScore;
                chainPred[it->Index] = bestSeed->Index;

                if (bestScore > bestChainScore) {
                    bestChainScore = bestScore;
                    bestChainEnd = it->Index;
                }
            } else {
                // PLEASE NOTE: these have already been done at creation time
                // scores[it->Index] = seedSize(*it);
                // chainPred[it->Index] = boost::none;
            }
        }

        sweepSet.insert(start, it);

        // remove all seeds from the sweepSet with end position less than the current row,
        // and ensure the colSet invariant is kept:
        //   that all columns greater than our current
        for (; toRemove != seeds.end() && (*toRemove).EndPositionV() < row; ++toRemove) {
            SDPColumn col((*toRemove).EndPositionH(), boost::make_optional(*toRemove));

            auto it = colSet.find(col);

            // update the column if it doesn't exist
            // or if its score is less than the fragment we're removing from consideration
            if (it == colSet.end() || zScore(*(it->Seed)) < zScore(*toRemove)) {
                // insert the updated column, get successor
                it = next(colSet.insert(col).first);

                // keep removing columns long as the scores are less than
                while (it != colSet.end() && zScore(*(it->Seed)) < zScore(*toRemove)) {
                    it = colSet.erase(it);
                }
            }

            sweepSet.erase(*toRemove);
        }
    }

    // seeds need to be sorted by Index to ... index into it properly
    sort(seeds.begin(), seeds.end(), IndexCompare);

    std::vector<PacBio::Align::Seed> chain;
    while (bestChainEnd) {
        chain.push_back(seeds[*bestChainEnd]);
        bestChainEnd = chainPred[*bestChainEnd];
    }

    // Seeds were added back-to-front, so reverse the current order in place
    std::reverse(chain.begin(), chain.end());

    return chain;
}

}  // namespace CCS
}  // namespace PacBio
