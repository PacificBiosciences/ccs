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

#include <seqan/seeds.h>
#include <seqan/sequence.h>

#include <pacbio/ccs/ChainSeeds.h>

namespace PacBio {
namespace CCS {
namespace {

// H is query, V is reference

long Diagonal(const seqan::Seed<seqan::Simple>& seed)
{
    size_t v = beginPositionV(seed);
    size_t h = beginPositionH(seed);

    if (v > h)
        return -static_cast<long>(v - h);

    return static_cast<long>(h - v);
}

bool HVCompare(const seqan::Seed<seqan::Simple>& a, const seqan::Seed<seqan::Simple>& b)
{
    unsigned leftA = beginPositionH(a),
             leftB = beginPositionH(b);

    if (leftA < leftB)
        return true;

    if (leftA == leftB)
        return endPositionV(a) < endPositionV(b);

    return false;
}

bool VHCompare(const seqan::Seed<seqan::Simple>& a, const seqan::Seed<seqan::Simple>& b)
{
    unsigned leftA = beginPositionV(a),
             leftB = beginPositionV(b);

    if (leftA < leftB)
        return true;

    if (leftA == leftB)
        return endPositionH(a) < endPositionH(b);

    return false;
}

bool DiagonalCompare(const seqan::Seed<seqan::Simple>& a, const seqan::Seed<seqan::Simple>& b)
{
    int diagA = Diagonal(a),
        diagB = Diagonal(b);

    if (diagA == diagB)
        return beginPositionH(a) < beginPositionH(b);

    return diagA < diagB;
}

long LinkScore(const seqan::Seed<seqan::Simple>& a, const seqan::Seed<seqan::Simple>& b, const int matchReward)
{
    using namespace std;

    long aH = static_cast<long>(beginPositionH(a));
    long aV = static_cast<long>(beginPositionV(a));
    long bH = static_cast<long>(beginPositionH(b));
    long bV = static_cast<long>(beginPositionV(b));
    long k  = min(static_cast<long>(seedSize(a)), static_cast<long>(seedSize(b)));

    long diagA = Diagonal(a);
    long diagB = Diagonal(b);
    long fwd = min(aH - bH, aV - bV);
    long indels = abs(diagA - diagB);
    long matches = k - max(0L, k - fwd);
    long mismatches = fwd - matches;

    return matchReward * matches - indels - mismatches;
}

struct SDPHit : public seqan::Seed<seqan::Simple>
{
    size_t Index;

    SDPHit(const seqan::Seed<seqan::Simple>& seed, const size_t index)
        : seqan::Seed<seqan::Simple>(seed)
        , Index(index)
    {}

    bool operator<(const SDPHit& other) const
    {
        return DiagonalCompare(*this, other);
    }
};

bool IndexCompare(const SDPHit& a, const SDPHit& b)
{
    return a.Index < b.Index;
}

struct SDPColumn
{
    boost::optional<SDPHit> Seed;
    size_t Column;

    SDPColumn(size_t column, const boost::optional<SDPHit> seed = boost::none)
        : Seed(seed)
        , Column(column)
    {}

    bool operator<(const SDPColumn& other) const
    {
        return Column < other.Column;
    }
};

std::vector<boost::optional<SDPHit>>
ComputeVisibilityLeft(const std::vector<SDPHit>& seeds,
                      std::set<SDPHit>& sweepSet)
{
    using namespace seqan;
    using namespace std;

    vector<boost::optional<SDPHit>> visible(seeds.size(), boost::none);
    auto toRemove = seeds.begin();

    for (auto it = seeds.begin(); it != seeds.end(); )
    {
        const size_t col = beginPositionH(*it);
        const auto start = it;

        for (; it != seeds.end() && col == beginPositionH(*it); ++it)
        {
            auto succ = sweepSet.upper_bound(*it);

            if (succ != sweepSet.end())
            {
                visible[it->Index] = *succ;
            }
        }

        sweepSet.insert(start, it);

        for (; toRemove != seeds.end() && endPositionH(*toRemove) < col; ++toRemove)
        {
            sweepSet.erase(*toRemove);
        }
    }

    // clear the sweepSet after use
    sweepSet.clear();

    return visible;
}

} // anonymous namespace


void
ChainSeeds(seqan::String<seqan::Seed<seqan::Simple>>& chain,
           const seqan::SeedSet<seqan::Seed<seqan::Simple>>& seedSet,
           const int matchReward)
{
    using namespace seqan;
    using namespace std;

    vector<SDPHit> seeds;
    vector<long> scores(length(seedSet), 0L);

    // initialize our "SDPHits" vector (has fixed index)
    {
        size_t i = 0;

        for (auto it = begin(seedSet); it != end(seedSet); ++it, ++i)
        {
            seeds.push_back(SDPHit(*it, i));
            scores[i] = seedSize(*it);
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

    auto zScore = [&scores] (const SDPHit& seed)
                  { return scores[seed.Index] + beginPositionH(seed) + beginPositionV(seed); };

    for (auto it = seeds.begin(); it != seeds.end(); )
    {
        const size_t row = beginPositionV(*it);
        const auto start = it;

        for (; it != seeds.end() && row == beginPositionV(*it); ++it)
        {
            long bestScore = -numeric_limits<long>::max();
            boost::optional<SDPHit> bestSeed = boost::none;

            // find the previous column and best fragment from it
            {
                SDPColumn col(beginPositionH(*it));
                auto pred = colSet.lower_bound(col);

                if (pred != colSet.begin())
                {
                    pred = prev(pred);
                    long s = scores[pred->Seed->Index] + LinkScore(*it, *(pred->Seed), matchReward);

                    if (s > bestScore)
                    {
                        bestScore = s;
                        bestSeed  = pred->Seed;
                    }
                }
            }

            // search visible fragments (above)
            {
                auto visa = sweepSet.lower_bound(*it);

                if (visa != sweepSet.begin())
                {
                    visa = prev(visa);
                    long s = scores[visa->Index] + LinkScore(*it, *visa, matchReward);

                    if (s > bestScore)
                    {
                        bestScore = s;
                        bestSeed  = *visa;
                    }
                }
            }

            // search visible fragments (left)
            if (auto visl = visible[it->Index])
            {
                long s = scores[visl->Index] + LinkScore(*it, *visl, matchReward);

                if (s > bestScore)
                {
                    bestScore = s;
                    bestSeed  = visl;
                }
            }

            if (bestSeed && bestScore > 0)
            {
                scores[it->Index] = bestScore;
                chainPred[it->Index] = bestSeed->Index;

                if (bestScore > bestChainScore)
                {
                    bestChainScore = bestScore;
                    bestChainEnd   = it->Index;
                }
            }
            else
            {
                // PLEASE NOTE: these have already been done at creation time
                // scores[it->Index] = seedSize(*it);
                // chainPred[it->Index] = boost::none;
            }
        }

        sweepSet.insert(start, it);

        // remove all seeds from the sweepSet with end position less than the current row,
        // and ensure the colSet invariant is kept:
        //   that all columns greater than our current
        for (; toRemove != seeds.end() && endPositionV(*toRemove) < row; ++toRemove)
        {
            SDPColumn col(endPositionH(*toRemove), boost::make_optional(*toRemove));

            auto it = colSet.find(col);

            // update the column if it doesn't exist
            // or if its score is less than the fragment we're removing from consideration
            if (it == colSet.end() || zScore(*(it->Seed)) < zScore(*toRemove))
            {
                // insert the updated column, get successor
                it = next(colSet.insert(col).first);

                // keep removing columns long as the scores are less than
                while (it != colSet.end() && zScore(*(it->Seed)) < zScore(*toRemove))
                {
                    it = colSet.erase(it);
                }
            }

            sweepSet.erase(*toRemove);
        }
    }

    // seeds need to be sorted by Index to ... index into it properly
    sort(seeds.begin(), seeds.end(), IndexCompare);
    clear(chain);

    while (bestChainEnd)
    {
        appendValue(chain, seeds[*bestChainEnd]);
        bestChainEnd = chainPred[*bestChainEnd];
    }

    // reverse the chain order in place
    reverse(chain);
}

} // namespace CCS
} // namespace PacBio
