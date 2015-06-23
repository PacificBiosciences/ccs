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

// Author: Lance Hepler, Brett Bowman

#pragma once

#include <algorithm>
#include <set>
#include <queue>

#include <boost/optional.hpp>
#include <boost/optional/optional_io.hpp>

#include <seqan/seeds.h>
#include <seqan/sequence.h>

#include <pbsparse/ChainSeedsConfig.h>

namespace PacBio {
namespace SparseAlignment {

// H is query, V is reference

namespace {

/// Calculate number of bases between a seed and the diagonal
///  axis of the matrix it's in as extended outward from the 
///  upper-left-hand corner of the matrix toward the lower right.
///  High numbers are closer to the upper-right, negative numbers
///  closer to the lower-left.
///
/// \param  seed  The seed whose diagonal is desired.
///
/// \return  long  The distance from the seed's upper-left corner and 
///                the diagonal axis of the matrix.
long Diagonal(const seqan::Seed<seqan::Simple>& seed)
{
    size_t v = beginPositionV(seed);
    size_t h = beginPositionH(seed);

    if (v > h)
        return -static_cast<long>(v - h);

    return static_cast<long>(h - v);
}

/// Compare seeds for sorting, first in the horizontal (query)
///  dimension, then in the vertical (reference) dimension.
///
/// \param  lhs  The left-hand side seed to compare
/// \param  rhs  The right-hand side seed to compare
///
/// \return  bool  Whether the left seed preceeds the right seed
bool HVCompare(const seqan::Seed<seqan::Simple>& lhs, 
               const seqan::Seed<seqan::Simple>& rhs)
{
    unsigned leftH  = beginPositionH(lhs),
             rightH = beginPositionH(rhs);

    if (leftH < rightH)
        return true;

    if (leftH == rightH)
        return endPositionV(lhs) < endPositionV(rhs);

    return false;
}

/// Compare seeds for sorting, first in the vertical (reference)
///  dimension, then in the horizontal (query) dimension.
///
/// \param  lhs  The left-hand side seed to compare
/// \param  rhs  The right-hand side seed to compare
///
/// \return  bool  Whether the left seed preceeds the right seed
bool VHCompare(const seqan::Seed<seqan::Simple>& lhs, 
               const seqan::Seed<seqan::Simple>& rhs)
{
    unsigned leftV  = beginPositionV(lhs),
             rightV = beginPositionV(rhs);

    if (leftV < rightV)
        return true;

    if (leftV == rightV)
        return endPositionH(lhs) < endPositionH(rhs);

    return false;
}

/// Compare seeds for sorting, by whether one seed is higher or
/// lower than another in the sparse matrix of their alignment, 
/// according to their diagonals.  Seeds near the upper-right
/// corner are said to preceed seeds closer to the lower-left.
///
/// \param  lhs  The left-hand seed to compare
/// \param  rhs  The right-hand seed to compare
///
/// \return  bool  Whether the left seed preceeds the right seed
bool DiagonalCompare(const seqan::Seed<seqan::Simple>& lhs, 
                     const seqan::Seed<seqan::Simple>& rhs)
{
    int leftD  = beginDiagonal(lhs),
        rightD = beginDiagonal(rhs);

    if (leftD == rightD)
        return beginPositionH(lhs) < beginPositionH(rhs);

    return leftD < rightD;
}

/// Score the possible linkage of two seeds based on 3 criteria:
///     (A) The number of bases in the shortest seed
///     (B) The number of bases between the two seeds
///     (C) The size of the difference between their diagonals
/// each with it's own weight(s).
///
/// \param  lhs  The left-hand side seed to link
/// \param  rhs  The right-hand side seed to link
/// \param  matchScore  The reward for each base within the seeds.
///                     (Value should be > 0).
/// \param  mismatchScore  The penalty for each base between the seeds.
///                        (Value should be <= 0)
/// \param  insertionScore  The penalty for each difference in between the
///                         diagonals, if Left > Right. (Value should be < 0)
/// \param  deletionScore  The penalty for each difference in between the
///                        diagonals, if Right > Left. (Value should be < 0)
///
/// \return  long  The score associated with the linkage
long LinkScore(const seqan::Seed<seqan::Simple>& lhs, 
               const seqan::Seed<seqan::Simple>& rhs, 
               const ChainSeedsConfig config)
{
    using namespace std;

    long lH  = static_cast<long>(beginPositionH(lhs));
    long lV  = static_cast<long>(beginPositionV(lhs));
    long rH  = static_cast<long>(beginPositionH(rhs));
    long rV  = static_cast<long>(beginPositionV(rhs));
    long k   = min(static_cast<long>(seedSize(lhs)), 
                   static_cast<long>(seedSize(rhs)));
    long fwd = min(lH - rH, lV - rV);

    // matchReward = # of anchor bases * matchScore;
    long matches = k - max(0l, k - fwd);
    long matchReward = matches * config.matchScore;
    
    // nonMatchPenalty = # of non-anchor, on-diagonal bases * nonMatchPenalty
    long nonMatches = fwd - matches;
    long nonMatchScorePenalty = nonMatches * config.nonMatchPenalty;
    
    // indelPenalty = difference in the seed diagonals * indelScore
    long diagL = Diagonal(lhs);
    long diagR = Diagonal(rhs);
    long drift = diagL - diagR;
    long indelScorePenalty = 0;
    if (drift > 0)
        indelScorePenalty = drift * config.insertionPenalty;
    else if (drift < 0)
        indelScorePenalty = -drift * config.deletionPenalty;

    return matchReward + indelScorePenalty + nonMatchScorePenalty;
}

///
/// A Sparse Dynamic Programming hit.  A wrapper around SeqAn's Seed class
/// with an additional field for storing it's index in the original seed set.
///
/// TODO (bbowman): Can I replace this struct with just Seeds by using the
///                 build in SeedScore field and setter/getters?
///
struct SDPHit : public seqan::Seed<seqan::Simple>
{
    size_t Index;

    SDPHit(const seqan::Seed<seqan::Simple>& seed, const size_t index)
        : seqan::Seed<seqan::Simple>(seed)
        , Index(index)
    {}

    ///
    ///
    /// \return  bool  Whether this s
    bool operator<(const SDPHit& other) const
    {
        return DiagonalCompare(*this, other);
    }
};

/// Compare two seeds (SDPHits) according to their original indices in
/// their source seedSet
///
/// \param  lhs  The left-hand side seed to compare
/// \param  rhs  The right-hand side seed to compare
///
/// \return  bool  Whether the original index of LHS is less than RHS
bool IndexCompare(const SDPHit& lhs, const SDPHit& rhs)
{
    return lhs.Index < rhs.Index;
}

///
/// A column in the Sparse Dynamic Programming matrix.  A wrapper around the above
///  SDPHit struct to give it a column field and a comparison operator that uses it.
/// 
/// TODO (bbowman): Can I replace this struct with just SDPHits or just Seeds, since
///                 we get the Column from SeqAn function endPositionH() function
///                 anyway?
///
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

/// 
/// 
/// \param  seeds  A vector of the SDPHits, sorted by their position in the
///                horizontal (query)
/// \param  sweepSet  
///
/// \return  vector<optional<SDPHits>>  For each seed in the input vector, the
///             first index that is visible to it's left, if any.
std::vector<boost::optional<SDPHit>>
ComputeVisibilityLeft(const std::vector<SDPHit>& seeds,
                      std::set<SDPHit>& sweepSet)
{
    using namespace seqan;
    using namespace std;
    using boost::optional;

    vector<optional<SDPHit>> visible(seeds.size(), boost::none);  // Output

    auto toRemove = seeds.begin();
    for (auto it = seeds.begin(); it != seeds.end(); )
    {
        const size_t col = beginPositionH(*it);
        const auto start = it;

        // Advance the iterator to the end of the current column in our
        //  column-sorted vector of seeds
        for (; it != seeds.end() && col == beginPositionH(*it); ++it)
        {
            // For each seed, record in the output vector, record the first
            //  seed after it in the sweepSet (if any).  Since the sweepSet
            //  only contains seeds from previous columns and is sorted by 
            //  their diagonals, this means the seeds found in this way will
            //  all (A) Start to the left and (B) Start on a higher diagonal
            auto succ = sweepSet.upper_bound(*it);
            if (succ != sweepSet.end())
            {
                visible[it->Index] = *succ;
            }
        }

        // Add all seeds to the sweepSet that start in the current column
        sweepSet.insert(start, it);

        // Remove all seeds from the sweepSet that end before the current column
        for (; toRemove != seeds.end() && endPositionH(*toRemove) < col; ++toRemove)
        {
            sweepSet.erase(*toRemove);
        }
    }

    // clear the sweepSet after use
    sweepSet.clear();

    return visible;
}

///
/// A possible chain of SDP seeds.  A simple struct for wrapping the
///  three pieces of information we need to filter and possibly 
///  reconstruct a chain that we've found:
///     (A) The reference sequence where the chain was found
///     (B) The terminal seed in the chain
///     (C) The chain's score
///
struct ChainHit
{
    size_t seedSetIdx;
    size_t endIndex;
    long score;
};

//
// Functor for sorting ChainHits by score
//
class ChainHitCompare
{
public:
    bool operator()(const ChainHit& lhs,
                    const ChainHit& rhs)
    {
        return lhs.score > rhs.score;
    }
};

/// Though we expect to receive the Seeds we'll be chaining in a tree-like 
///  SeedSet, we need them and their scores in well-ordered vectors to
///  perform the actual chaining ourselves.  This function is a simple
///  helper for abstracting away that process.
///
/// \param  seedSet  The collection of seeds to be chained
/// \param  seeds    The vector of SDPHit objects that will actually be chained
/// \param  scores   The vector of scores to initialize with default per-seed scores
void InitializeSeedsAndScores(const seqan::SeedSet<seqan::Seed<seqan::Simple>>& seedSet,
                              std::vector<SDPHit>* seeds,
                              std::vector<long>* scores)
{
    size_t i = 0;
    for (auto it = begin(seedSet); it != end(seedSet); ++it, ++i)
    {
        seeds->push_back(SDPHit(*it, i));
        scores->at(i) = seedSize(*it);
    }
}


/// Search a SeedSet for the best numCandidates sets of locally-chainable 
/// seeds according to some scoring criteria.  Seed chains are scored based
/// on their length and penalized according to the distance between them and
/// how far apart their diagonals are.  Final scores for a chain must be above
/// some minScore threshold to be reported.
///
/// Roughly equivalent to BLASR's SDPAlign
///
/// TODO (bbowman): I shouldn't report partial and complete chains of the same seeds
/// TODO (bbowman): Figure out why my penalties need to be lower than BLASR's
///                 for similar results
///
/// \param  chains  A vector of SeedSets to store and return collections of 
///                 locally chained seeds in.
/// \param  seedSet  The SeedSet to search for chains in
/// \param  numCandidates  The maximum number of chains to report
/// \param  minScore  The minimum score to require from a reported chain
/// \param  match  Score to add to each matching base in each seed
/// \param  nonmatch  Penalty for each non-aligned base between seeds
///                   along whichever dimension they are closest together
/// \param  insertion  Penalty for each base along the diagonal that
///                    separates two seeds, if the 'upstream' seed is on 
///                    a higher diagonal than the 'downstream' seed.
/// \param  deletion  Penalty for each base along the diagonal that
///                   separates two seeds, if the 'upstream' seed is on 
///                   a lower diagonal than the 'downstream' seed.
void __attribute__((__unused__))
ChainSeedsImpl(std::priority_queue<ChainHit, std::vector<ChainHit>, ChainHitCompare>* chainHits,
                    std::vector<boost::optional<size_t>>* chainPred,
                    std::vector<SDPHit>* seeds,
                    std::vector<long>& scores,
                    const size_t seedSetIdx,
                    const ChainSeedsConfig& config)
{
    using namespace std;
    using namespace seqan;

    // compute visibility left, requires H-sorted seeds
    std::set<SDPHit> sweepSet;
    sort(seeds->begin(), seeds->end(), HVCompare);
    auto visible = ComputeVisibilityLeft(*seeds, sweepSet);

    // compute the visibility above, requires V-sorted seeds
    sort(seeds->begin(), seeds->end(), VHCompare);
    auto toRemove = seeds->begin();
    std::set<SDPColumn> colSet;

    auto zScore = [&scores] (const SDPHit& seed)
                  { return scores[seed.Index] + beginPositionH(seed) + beginPositionV(seed); };

    for (auto it = seeds->begin(); it != seeds->end(); )
    {
        const size_t row = beginPositionV(*it);
        const auto start = it;

        for (; it != seeds->end() && row == beginPositionV(*it); ++it)
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
                    long s = scores[pred->Seed->Index] + LinkScore(*it, *(pred->Seed), config);

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
                    long s = scores[visa->Index] + LinkScore(*it, *visa, config);

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
                long s = scores[visl->Index] + LinkScore(*it, *visl, config);

                if (s > bestScore)
                {
                    bestScore = s;
                    bestSeed  = visl;
                }
            }

            if (bestSeed && bestScore >= config.minScore)
            {
                scores[it->Index] = bestScore;
                chainPred->at(it->Index) = bestSeed->Index;

                if (chainHits->size() < config.numCandidates)
                {
                    chainHits->push( {seedSetIdx, it->Index, bestScore} );
                }
                else if (bestScore > chainHits->top().score)
                {
                    chainHits->pop();
                    chainHits->push( {seedSetIdx, it->Index, bestScore} );
                }
            }
            else if (scores[it->Index] >= config.minScore)
            {
                // PLEASE NOTE: these have already been done at creation time
                // scores[it->Index] = seedSize(*it);
                // chainPred[it->Index] = boost::none;
                //

                if (chainHits->size() < config.numCandidates)
                {
                    chainHits->push( {seedSetIdx, it->Index, bestScore} );
                }
                else if (bestScore > chainHits->top().score)
                {
                    chainHits->pop();
                    chainHits->push( {seedSetIdx, it->Index, bestScore} );
                }
            }
        }

        sweepSet.insert(start, it);

        // remove all seeds from the sweepSet with end position less than the current row,
        // and ensure the colSet invariant is kept:
        //   that all columns greater than our current
        for (; toRemove != seeds->end() && endPositionV(*toRemove) < row; ++toRemove)
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
    sort(seeds->begin(), seeds->end(), IndexCompare);
}

/// Search a SeedSet for the best numCandidates sets of locally-chainable 
/// seeds according to some scoring criteria.  Seed chains are scored based
/// on their length and penalized according to the distance between them and
/// how far apart their diagonals are.  Final scores for a chain must be above
/// some minScore threshold to be reported.
///
/// Roughly equivalent to BLASR's SDPAlign
///
/// \param  chains  A vector of SeedSets to store and return collections of 
///                 locally chained seeds in.
/// \param  seedSet  The SeedSet to search for chains in
/// \param  numCandidates  The maximum number of chains to report
/// \param  minScore  The minimum score to require from a reported chain
/// \param  match  Score to add to each matching base in each seed
/// \param  nonmatch  Penalty for each non-aligned base between seeds
///                   along whichever dimension they are closest together
/// \param  insertion  Penalty for each base along the diagonal that
///                    separates two seeds, if the 'upstream' seed is on 
///                    a higher diagonal than the 'downstream' seed.
/// \param  deletion  Penalty for each base along the diagonal that
///                   separates two seeds, if the 'upstream' seed is on 
///                   a lower diagonal than the 'downstream' seed.
void __attribute__((__unused__))
ChainSeeds(std::vector<seqan::String<seqan::Seed<seqan::Simple>>>* chains,
                const seqan::SeedSet<seqan::Seed<seqan::Simple>>& seedSet,
                const ChainSeedsConfig& config)
{
    using namespace seqan;
    using namespace std;

    // Initialize the work-horse vectors we will actually work with
    priority_queue<ChainHit, vector<ChainHit>, ChainHitCompare> chainHits;    
    vector<boost::optional<size_t>> chainPred(length(seedSet), boost::none);
    vector<SDPHit> seeds;
    vector<long> scores(length(seedSet), 0L);
    InitializeSeedsAndScores(seedSet, &seeds, &scores);

    // Call the main function
    ChainSeedsImpl(&chainHits, &chainPred, &seeds, scores, 0, config);

    // Empty and resize the result vector
    chains->clear();
    chains->resize(chainHits.size());

    // Pop our results from our queue and convert them into Seed Chains / Sets
    int i = chainHits.size() - 1;
    while (!chainHits.empty())
    {
        const auto hit = chainHits.top();
        //std::cout << "b(" << hit.reference << ", " << hit.endIndex << ", " << hit.score << ")" << std::endl;
        
        // While there are additional links in the chain, append them
        boost::optional<size_t> chainEnd = hit.endIndex;
        while (chainEnd)
        {
            appendValue(chains->at(i), seeds[*chainEnd]);
            chainEnd = chainPred[*chainEnd];
        }
        
        // We appended seeds back-to-front, so reverse the current order in place
        reverse(chains->at(i));

        chainHits.pop();
        --i;
    }
}


/// Search a SeedSet for the best numCandidates sets of locally-chainable 
/// seeds according to some scoring criteria.  Seed chains are scored based
/// on their length and penalized according to the distance between them and
/// how far apart their diagonals are.  Final scores for a chain must be above
/// some minScore threshold to be reported.
///
/// Roughly equivalent to BLASR's SDPAlign
///
/// \param  chains  A vector of SeedSets to store and return collections of 
///                 locally chained seeds in.
/// \param  seedSet  The SeedSet to search for chains in
/// \param  numCandidates  The maximum number of chains to report
/// \param  minScore  The minimum score to require from a reported chain
/// \param  match  Score to add to each matching base in each seed
/// \param  nonmatch  Penalty for each non-aligned base between seeds
///                   along whichever dimension they are closest together
/// \param  insertion  Penalty for each base along the diagonal that
///                    separates two seeds, if the 'upstream' seed is on 
///                    a higher diagonal than the 'downstream' seed.
/// \param  deletion  Penalty for each base along the diagonal that
///                   separates two seeds, if the 'upstream' seed is on 
///                   a lower diagonal than the 'downstream' seed.
void __attribute__((__unused__))
ChainSeeds(std::vector<seqan::SeedSet<seqan::Seed<seqan::Simple>>>* chains,
                const seqan::SeedSet<seqan::Seed<seqan::Simple>>& seedSet,
                const ChainSeedsConfig& config)
{
    using namespace seqan;
    using namespace std;

    // Initialize the work-horse vectors we will actually work with
    priority_queue<ChainHit, vector<ChainHit>, ChainHitCompare> chainHits;    
    vector<boost::optional<size_t>> chainPred(length(seedSet), boost::none);
    vector<SDPHit> seeds;
    vector<long> scores(length(seedSet), 0L);
    InitializeSeedsAndScores(seedSet, &seeds, &scores);

    // Call the main function
    ChainSeedsImpl(&chainHits, &chainPred, &seeds, scores, 0, config);

    // Empty and resize the result vector
    chains->clear();
    chains->resize(chainHits.size());

    // Pop our results from our queue and convert them into Seed Chains / Sets
    int i = chainHits.size() - 1;
    while (!chainHits.empty())
    {
        const auto hit = chainHits.top();
        
        // While there are additional links in the chain, append them
        boost::optional<size_t> chainEnd = hit.endIndex;
        while (chainEnd)
        {
            Seed<Simple> seed = seeds[*chainEnd];
            addSeed(chains->at(i), seed, Single());
            chainEnd = chainPred[*chainEnd];
        }
        
        chainHits.pop();
        --i;
    }
}
    
/// Search a SeedSet for the best numCandidates sets of locally-chainable 
/// seeds according to some scoring criteria.  Seed chains are scored based
/// on their length and penalized according to the distance between them and
/// how far apart their diagonals are.  Final scores for a chain must be above
/// some minScore threshold to be reported.
///
/// Roughly equivalent to BLASR's SDPAlign
///
/// \param  chains  A vector of SeedSets to store and return collections of 
///                 locally chained seeds in.
/// \param  seedSet  The SeedSet to search for chains in
/// \param  numCandidates  The maximum number of chains to report
/// \param  minScore  The minimum score to require from a reported chain
/// \param  match  Score to add to each matching base in each seed
/// \param  nonmatch  Penalty for each non-aligned base between seeds
///                   along whichever dimension they are closest together
/// \param  insertion  Penalty for each base along the diagonal that
///                    separates two seeds, if the 'upstream' seed is on 
///                    a higher diagonal than the 'downstream' seed.
/// \param  deletion  Penalty for each base along the diagonal that
///                   separates two seeds, if the 'upstream' seed is on 
///                   a lower diagonal than the 'downstream' seed.
void ChainSeeds(std::vector<std::pair<size_t, seqan::SeedSet<seqan::Seed<seqan::Simple>>>>* chains,
                const std::map<size_t, seqan::SeedSet<seqan::Seed<seqan::Simple>>> seedSets,
                const ChainSeedsConfig config)
{
    using namespace seqan;
    using namespace std;

    // The queue will accumulat results across SeedSets
    priority_queue<ChainHit, vector<ChainHit>, ChainHitCompare> chainHits;

    // Our vectors of seeds and chain-links need to persist for eventual
    //  use reconstructing our chains, so we initialize them as 
    //  vectors-of-vectors here, 1 per seedSet we will analyze
    size_t numSeedSets = seedSets.size();
    vector<vector<boost::optional<size_t>>> chainPred(numSeedSets);
    vector<vector<SDPHit>> seeds(numSeedSets);

    // We also need to record which seedSet came from which reference
    vector<size_t> references(numSeedSets);

    // Iterate over the multiple SeedSets once to search for chains
    int i = 0;
    for (auto it = seedSets.begin(); it != seedSets.end(); ++it, ++i)
    {
        // Extract the contents of our K/V pair
        references[i] = it->first;
        const auto& seedSet = it->second;

        // Initialize the work-horse vectors we will actually work with
        chainPred[i] = vector<boost::optional<size_t>>(length(seedSet), boost::none);
        vector<long> scores(length(seedSet), 0L);
        InitializeSeedsAndScores(seedSet, &seeds[i], &scores);

        // Call the main function on the current seedSet
        ChainSeedsImpl(&chainHits, &chainPred[i], &seeds[i], scores, i, config);
    }

    // Empty and resize the result vector
    chains->clear();
    chains->resize(chainHits.size());

    // Pop our results from our queue and convert them into Seed Chains / Sets
    int j = chainHits.size() - 1;
    while (!chainHits.empty())
    {
        // Take our next hit and extract the relevant indices ...
        const auto hit = chainHits.top();
        chains->at(j).first = references[hit.seedSetIdx];

        // While there are additional links in the chain, append them
        boost::optional<size_t> chainEnd = hit.endIndex;
        while (chainEnd)
        {
            Seed<Simple> seed = seeds[hit.seedSetIdx][*chainEnd];
            addSeed(chains->at(j).second, seed, Single());
            chainEnd = chainPred[hit.seedSetIdx][*chainEnd];
        }
        
        chainHits.pop();
        --j;
    }
}

} // anonymous namespace

}}  // PacBio::SparseAlignment
