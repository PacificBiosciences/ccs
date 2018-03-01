// Author: Brett Bowman

#pragma once

namespace PacBio {
namespace Align {

// Default values for struct fields
static const size_t NUM_CANDIDATES = 10;
static const long MIN_SCORE = 18;
static const int MATCH_SCORE = 5;
static const int NON_MATCH_PENALTY = 0;
static const int INSERTION_PENALTY = -4;
static const int DELETION_PENALTY = -8;
static const int MAX_SEED_GAP = 200;

/// A simple struct for represting a complete specialization for an
///  implementation of the Baker-Giancarlo SDP algorithm.  Templating
///  around this struct allows us to greatly reduce the number of
///  parameters that need to be thrown around.
///
struct ChainSeedsConfig
{
    // Default constructor
    ChainSeedsConfig(const size_t numCandidatesArg, const long minScoreArg, const int matchScoreArg,
                     const int nonMatchPenaltyArg, const int insertionPenaltyArg,
                     const int deletionPenaltyArg, const int maxSeedGapArg)
        : numCandidates(numCandidatesArg)
        , minScore(minScoreArg)
        , matchScore(matchScoreArg)
        , nonMatchPenalty(nonMatchPenaltyArg)
        , insertionPenalty(insertionPenaltyArg)
        , deletionPenalty(deletionPenaltyArg)
        , maxSeedGap(maxSeedGapArg)
    {
    }

    // Abbreviated constructors
    ChainSeedsConfig(const size_t numCandidatesArg)
        : numCandidates(numCandidatesArg)
        , minScore(MIN_SCORE)
        , matchScore(MATCH_SCORE)
        , nonMatchPenalty(NON_MATCH_PENALTY)
        , insertionPenalty(INSERTION_PENALTY)
        , deletionPenalty(DELETION_PENALTY)
        , maxSeedGap(MAX_SEED_GAP)
    {
    }

    ChainSeedsConfig()
        : numCandidates(NUM_CANDIDATES)
        , minScore(MIN_SCORE)
        , matchScore(MATCH_SCORE)
        , nonMatchPenalty(NON_MATCH_PENALTY)
        , insertionPenalty(INSERTION_PENALTY)
        , deletionPenalty(DELETION_PENALTY)
        , maxSeedGap(MAX_SEED_GAP)
    {
    }

    size_t numCandidates;
    long minScore;
    int matchScore;
    int nonMatchPenalty;
    int insertionPenalty;
    int deletionPenalty;
    int maxSeedGap;
};
}
}  // PacBio::Align
