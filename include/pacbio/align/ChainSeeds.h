// Authors: Lance Hepler, Brett Bowman

#pragma once

#include <algorithm>
#include <queue>
#include <set>

#include <boost/optional.hpp>
#include <boost/optional/optional_io.hpp>

#include <seqan/seeds.h>
#include <seqan/sequence.h>

#include <pbcopper/align/Seeds.h>

#include <pacbio/align/ChainSeedsConfig.h>

namespace PacBio {
namespace Align {

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
///
inline long Diagonal(const Seed& seed)
{
    const auto v = seed.BeginPositionV();
    const auto h = seed.BeginPositionH();

    if (v > h) return -static_cast<long>(v - h);

    return static_cast<long>(h - v);
}

/// Compare seeds for sorting, first in the horizontal (query)
///  dimension, then in the vertical (reference) dimension.
///
/// \param  lhs  The left-hand side seed to compare
/// \param  rhs  The right-hand side seed to compare
///
/// \return  bool  Whether the left seed preceeds the right seed
///
inline bool HVCompare(const Seed& lhs, const Seed& rhs)
{
    const auto leftH = lhs.BeginPositionH();
    const auto rightH = rhs.BeginPositionH();

    if (leftH < rightH) return true;

    if (leftH == rightH) return lhs.EndPositionV() < rhs.EndPositionV();

    return false;
}

/// Compare seeds for sorting, first in the vertical (reference)
///  dimension, then in the horizontal (query) dimension.
///
/// \param  lhs  The left-hand side seed to compare
/// \param  rhs  The right-hand side seed to compare
///
/// \return  bool  Whether the left seed preceeds the right seed
///
inline bool VHCompare(const Seed& lhs, const Seed& rhs)
{
    const auto leftV = lhs.BeginPositionV();
    const auto rightV = rhs.BeginPositionV();

    if (leftV < rightV) return true;

    if (leftV == rightV) return lhs.EndPositionH() < rhs.EndPositionH();

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
///
inline bool DiagonalCompare(const Seed& lhs, const Seed& rhs)
{
    const auto leftD = lhs.BeginDiagonal();
    const auto rightD = rhs.BeginDiagonal();

    if (leftD == rightD) return lhs.BeginPositionH() < rhs.BeginPositionH();

    return leftD < rightD;
}

/// Score the possible linkage of two seeds based on 3 criteria:
///     (A) The number of bases in the shortest seed
///     (B) The number of bases between the two seeds
///     (C) The size of the difference between their diagonals
/// each with it's own weight(s).
///
/// \param  lhs     The left-hand side seed to link
/// \param  rhs     The right-hand side seed to link
/// \param  config  Provides scoring values to use when chaining
///
/// \return  long  The score associated with the linkage
///
inline long LinkScore(const Seed& lhs, const Seed& rhs, const ChainSeedsConfig config)
{
    const auto lH = static_cast<long>(lhs.BeginPositionH());
    const auto lV = static_cast<long>(lhs.BeginPositionV());
    const auto rH = static_cast<long>(rhs.BeginPositionH());
    const auto rV = static_cast<long>(rhs.BeginPositionV());
    const auto k = std::min(static_cast<long>(lhs.Size()), static_cast<long>(rhs.Size()));
    const auto fwd = std::min(lH - rH, lV - rV);

    // matchReward = # of anchor bases * matchScore;
    const auto matches = k - std::max(0l, k - fwd);
    const auto matchReward = matches * config.matchScore;

    // nonMatchPenalty = # of non-anchor, on-diagonal bases * nonMatchPenalty
    const auto nonMatches = fwd - matches;
    const auto nonMatchScorePenalty = nonMatches * config.nonMatchPenalty;

    // Ignore any linkage over a certain size, no matter the score
    if (nonMatches > config.maxSeedGap) return -1;

    // indelPenalty = difference in the seed diagonals * indelScore
    const auto diagL = lhs.Diagonal();
    const auto diagR = rhs.Diagonal();
    const auto drift = diagL - diagR;
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
struct SDPHit : public Seed
{
    size_t Index;

    SDPHit(const Seed& seed, const size_t index) : Seed(seed), Index(index) {}

    ///
    ///
    /// \return  bool  Whether this s
    bool operator<(const SDPHit& other) const { return DiagonalCompare(*this, other); }
};

/// Compare two seeds (SDPHits) according to their original indices in
/// their source seedSet
///
/// \param  lhs  The left-hand side seed to compare
/// \param  rhs  The right-hand side seed to compare
///
/// \return  bool  Whether the original index of LHS is less than RHS
///
inline bool IndexCompare(const SDPHit& lhs, const SDPHit& rhs) { return lhs.Index < rhs.Index; }

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
        : Seed(seed), Column(column)
    {
    }

    bool operator<(const SDPColumn& other) const { return Column < other.Column; }
};

///
///
/// \param  seeds  A vector of the SDPHits, sorted by their position in the
///                horizontal (query)
/// \param  sweepSet
///
/// \return  vector<optional<SDPHits>>  For each seed in the input vector, the
///             first index that is visible to it's left, if any.
inline std::vector<boost::optional<SDPHit>> ComputeVisibilityLeft(const std::vector<SDPHit>& seeds,
                                                                  std::set<SDPHit>& sweepSet)
{
    using namespace seqan;
    using namespace std;
    using boost::optional;

    vector<optional<SDPHit>> visible(seeds.size(), boost::none);  // Output

    auto toRemove = seeds.begin();
    for (auto it = seeds.begin(); it != seeds.end();) {
        const auto col = (*it).BeginPositionH();
        const auto start = it;

        // Advance the iterator to the end of the current column in our
        //  column-sorted vector of seeds
        for (; it != seeds.end() && col == (*it).BeginPositionH(); ++it) {
            // For each seed, record in the output vector, record the first
            //  seed after it in the sweepSet (if any).  Since the sweepSet
            //  only contains seeds from previous columns and is sorted by
            //  their diagonals, this means the seeds found in this way will
            //  all (A) Start to the left and (B) Start on a higher diagonal
            auto succ = sweepSet.upper_bound(*it);
            if (succ != sweepSet.end()) {
                visible[it->Index] = *succ;
            }
        }

        // Add all seeds to the sweepSet that start in the current column
        sweepSet.insert(start, it);

        // Remove all seeds from the sweepSet that end before the current column
        for (; toRemove != seeds.end() && (*toRemove).EndPositionH() < col; ++toRemove) {
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
    bool operator()(const ChainHit& lhs, const ChainHit& rhs) { return lhs.score > rhs.score; }
};

/// Though we expect to receive the Seeds we'll be chaining in a tree-like
///  Seeds collection, we need them and their scores in well-ordered vectors to
///  perform the actual chaining ourselves.  This function is a simple
///  helper for abstracting away that process.
///
/// \param  seedSet  The collection of seeds to be chained
/// \param  seeds    The vector of SDPHit objects that will actually be chained
/// \param  scores   The vector of scores to initialize with default per-seed scores
inline void InitializeSeedsAndScores(const Seeds& seedSet, std::vector<SDPHit>* seeds,
                                     std::vector<long>* scores)
{
    size_t i = 0;
    for (const auto& seed : seedSet) {
        seeds->push_back(SDPHit(seed, i));
        scores->at(i) = seed.Size();
        ++i;
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
/// TODO (dbarnett): This parameter documentation below is _way_ out of sync.
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
inline void __attribute__((__unused__))
ChainSeedsImpl(std::priority_queue<ChainHit, std::vector<ChainHit>, ChainHitCompare>* chainHits,
               std::vector<boost::optional<size_t>>* chainPred, std::vector<SDPHit>* seeds,
               std::vector<long>& scores, const size_t seedSetIdx, const ChainSeedsConfig& config)
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

    auto zScore = [&scores](const SDPHit& seed) {
        return scores[seed.Index] + seed.BeginPositionH() + seed.BeginPositionV();
    };

    for (auto it = seeds->begin(); it != seeds->end();) {
        const size_t row = (*it).BeginPositionV();
        const auto start = it;

        for (; it != seeds->end() && row == (*it).BeginPositionV(); ++it) {
            long bestScore = -numeric_limits<long>::max();
            boost::optional<SDPHit> bestSeed = boost::none;

            // find the previous column and best fragment from it
            {
                SDPColumn col((*it).BeginPositionH());
                auto pred = colSet.lower_bound(col);

                if (pred != colSet.begin()) {
                    pred = prev(pred);
                    long s = scores[pred->Seed->Index] + LinkScore(*it, *(pred->Seed), config);

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
                    long s = scores[visa->Index] + LinkScore(*it, *visa, config);

                    if (s > bestScore) {
                        bestScore = s;
                        bestSeed = *visa;
                    }
                }
            }

            // search visible fragments (left)
            if (auto visl = visible[it->Index]) {
                long s = scores[visl->Index] + LinkScore(*it, *visl, config);

                if (s > bestScore) {
                    bestScore = s;
                    bestSeed = visl;
                }
            }

            if (bestSeed && bestScore >= config.minScore) {
                scores[it->Index] = bestScore;
                chainPred->at(it->Index) = bestSeed->Index;

                if (chainHits->size() < config.numCandidates) {
                    chainHits->push({seedSetIdx, it->Index, bestScore});
                } else if (bestScore > chainHits->top().score) {
                    chainHits->pop();
                    chainHits->push({seedSetIdx, it->Index, bestScore});
                }
            } else if (scores[it->Index] >= config.minScore) {
                // PLEASE NOTE: these have already been done at creation time
                // scores[it->Index] = seedSize(*it);
                // chainPred[it->Index] = boost::none;
                //

                if (chainHits->size() < config.numCandidates) {
                    chainHits->push({seedSetIdx, it->Index, bestScore});
                } else if (bestScore > chainHits->top().score) {
                    chainHits->pop();
                    chainHits->push({seedSetIdx, it->Index, bestScore});
                }
            }
        }

        sweepSet.insert(start, it);

        // remove all seeds from the sweepSet with end position less than the current row,
        // and ensure the colSet invariant is kept:
        //   that all columns greater than our current
        for (; toRemove != seeds->end() && (*toRemove).EndPositionV() < row; ++toRemove) {
            SDPColumn col((*toRemove).EndPositionH(), boost::make_optional(*toRemove));

            auto it = colSet.find(col);

            // update the column if it doesn't exist
            // or if its score is less than the fragment we're removing from consideration
            if (it == colSet.end() || zScore(*(it->Seed)) < zScore(*toRemove)) {
                // insert the updated column, get successor
                it = std::next(colSet.insert(col).first);

                // keep removing columns long as the scores are less than
                while (it != colSet.end() && zScore(*(it->Seed)) < zScore(*toRemove)) {
                    it = colSet.erase(it);
                }
            }

            sweepSet.erase(*toRemove);
        }
    }

    // seeds need to be sorted by Index to ... index into it properly
    sort(seeds->begin(), seeds->end(), IndexCompare);
}

/// Search a Seed set for the best numCandidates sets of locally-chainable
/// seeds according to some scoring criteria.  Seed chains are scored based
/// on their length and penalized according to the distance between them and
/// how far apart their diagonals are.  Final scores for a chain must be above
/// some minScore threshold to be reported.
///
/// Roughly equivalent to BLASR's SDPAlign
///
/// \param  seedSet The Seeds set to search for chains in
/// \param  config  Provides scoring values to use when chaining
///
/// \return  A vector of Seed vectors containing locally chained seeds.
///
inline std::vector<std::vector<Seed>> ChainSeeds(const Seeds& seedSet,
                                                 const ChainSeedsConfig& config)
{
    using namespace seqan;
    using namespace std;

    // Initialize the work-horse vectors we will actually work with
    priority_queue<ChainHit, vector<ChainHit>, ChainHitCompare> chainHits;
    vector<boost::optional<size_t>> chainPred(seedSet.size(), boost::none);
    vector<SDPHit> seeds;
    vector<long> scores(seedSet.size(), 0L);
    InitializeSeedsAndScores(seedSet, &seeds, &scores);

    // Call the main function
    ChainSeedsImpl(&chainHits, &chainPred, &seeds, scores, 0, config);

    // Empty and resize the result vector
    std::vector<std::vector<Seed>> chains;
    chains.resize(chainHits.size());

    // Pop our results from our queue and convert them into Seed Chains / Sets
    int i = chainHits.size() - 1;
    while (!chainHits.empty()) {
        const auto hit = chainHits.top();
        //std::cout << "b(" << hit.reference << ", " << hit.endIndex << ", " << hit.score << ")" << std::endl;

        // While there are additional links in the chain, append them
        boost::optional<size_t> chainEnd = hit.endIndex;
        while (chainEnd) {
            chains.at(i).push_back(seeds[*chainEnd]);
            chainEnd = chainPred[*chainEnd];
        }

        // We appended seeds back-to-front, so reverse the current order in place
        std::reverse(chains.at(i).begin(), chains.at(i).end());

        chainHits.pop();
        --i;
    }

    return chains;
}

///// Search a Seed set for the best numCandidates sets of locally-chainable
///// seeds according to some scoring criteria.  Seed chains are scored based
///// on their length and penalized according to the distance between them and
///// how far apart their diagonals are.  Final scores for a chain must be above
///// some minScore threshold to be reported.
/////
///// Roughly equivalent to BLASR's SDPAlign
/////
///// \param  seedSet The SeedSet to search for chains in
///// \param  config  Provides scoring values to use when chaining
/////
///// \return A vector of Seed sets containing locally chained seeds.
/////
inline std::vector<Seeds> ChainedSeedSets(const Seeds& seedSet, const ChainSeedsConfig& config)
{
    using namespace seqan;
    using namespace std;

    // Initialize the work-horse vectors we will actually work with
    priority_queue<ChainHit, vector<ChainHit>, ChainHitCompare> chainHits;
    vector<boost::optional<size_t>> chainPred(seedSet.size(), boost::none);
    vector<SDPHit> seeds;
    vector<long> scores(seedSet.size(), 0L);
    InitializeSeedsAndScores(seedSet, &seeds, &scores);

    // Call the main function
    ChainSeedsImpl(&chainHits, &chainPred, &seeds, scores, 0, config);

    // Empty and resize the result vector
    std::vector<Seeds> chains;
    chains.resize(chainHits.size());

    // Pop our results from our queue and convert them into Seed Chains / Sets
    int i = chainHits.size() - 1;
    while (!chainHits.empty()) {
        const auto hit = chainHits.top();

        // While there are additional links in the chain, append them
        boost::optional<size_t> chainEnd = hit.endIndex;
        while (chainEnd) {
            auto seed = seeds[*chainEnd];
            chains.at(i).AddSeed(seed);
            chainEnd = chainPred[*chainEnd];
        }

        chainHits.pop();
        --i;
    }
    return chains;
}

/// Search a Seed set for the best numCandidates sets of locally-chainable
/// seeds according to some scoring criteria.  Seed chains are scored based
/// on their length and penalized according to the distance between them and
/// how far apart their diagonals are.  Final scores for a chain must be above
/// some minScore threshold to be reported.
///
/// Roughly equivalent to BLASR's SDPAlign
///
/// \param  seedSet The SeedSet to search for chains in
/// \param  config  Provides scoring values to use when chaining
///
/// \return  A vector of SeedSets containing locally chained seeds.
///
inline std::vector<std::pair<size_t, Seeds>> ChainSeeds(const std::map<size_t, Seeds> seedSets,
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
    for (auto it = seedSets.begin(); it != seedSets.end(); ++it, ++i) {
        // Extract the contents of our K/V pair
        references[i] = it->first;
        const auto& seedSet = it->second;

        // Initialize the work-horse vectors we will actually work with
        chainPred[i] = vector<boost::optional<size_t>>(seedSet.size(), boost::none);
        vector<long> scores(seedSet.size(), 0L);
        InitializeSeedsAndScores(seedSet, &seeds[i], &scores);

        // Call the main function on the current seedSet
        ChainSeedsImpl(&chainHits, &chainPred[i], &seeds[i], scores, i, config);
    }

    // Empty and resize the result vector
    std::vector<std::pair<size_t, Seeds>> chains;
    chains.resize(chainHits.size());

    // Pop our results from our queue and convert them into Seed Chains / Sets
    int j = chainHits.size() - 1;
    while (!chainHits.empty()) {
        // Take our next hit and extract the relevant indices ...
        const auto hit = chainHits.top();
        chains.at(j).first = references[hit.seedSetIdx];

        // While there are additional links in the chain, append them
        boost::optional<size_t> chainEnd = hit.endIndex;
        while (chainEnd) {
            auto seed = seeds[hit.seedSetIdx][*chainEnd];
            chains.at(j).second.AddSeed(seed);
            chainEnd = chainPred[hit.seedSetIdx][*chainEnd];
        }

        chainHits.pop();
        --j;
    }

    return chains;
}

}  // anonymous namespace
}
}  // PacBio::Align
