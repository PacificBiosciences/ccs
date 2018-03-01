// Author: Brett Bowman

#pragma once

#include <seqan/index.h>
#include <seqan/seeds.h>

#include "SeedScorer.h"

namespace PacBio {
namespace Align {

template <typename TConfig>
class ChainScorer
{

public:  // Template types and values
    typedef typename TConfig::IndexType TIndexSpec;
    typedef typename TConfig::ShapeType TShapeSpec;
    size_t Size = TConfig::Size;

public:  // Derived types and values
    typedef seqan::Dna TAlphabet;
    typedef seqan::String<TAlphabet> TString;
    typedef seqan::Index<seqan::StringSet<TString>, TIndexSpec> TIndex;

    // Equivalent to BLASR's "T_Tuple"
    typedef seqan::Shape<TAlphabet, TShapeSpec> TShape;
    // Equivalent to BLASR's "MatchPos"
    typedef seqan::Seed<seqan::Simple> TSeed;
    // Equivalent to BLASR's "LongestIncreasingSubsequence"
    typedef seqan::String<TSeed> TSeedChain;
    typedef std::pair<size_t, TSeedChain> THit;

public:  // structors
    // Default constructor
    ChainScorer(const TIndex& index, const size_t kmerSize)
        : index_(index)
        , shape_(seqan::indexShape(index))
        , kmerSize_(kmerSize)
        , seedScorer_(index, kmerSize)
    {
        using namespace seqan;

        // Store the size of the reference as a float for calculating frequencies
        referenceSize_ = static_cast<float>(length(getFibre(index_, QGram_RawText()))) / 2.0f;
    }

    // Move constructor
    ChainScorer(ChainScorer&& src) = delete;
    // Copy constructor
    ChainScorer(const ChainScorer& src) = delete;
    // Move assignment operator
    ChainScorer& operator=(ChainScorer&& rhs) = delete;
    // Copy assignment operator
    ChainScorer& operator=(const ChainScorer& rhs) = delete;
    // Destructor
    ~ChainScorer() = default;

public:  // non-modifying methods
    /// Score a given seed-set found that matches the reference index, returning
    /// a score that approximates it's log-likehood.
    ///
    /// Equivalent to BLASR's ComputeLISPValue function.
    ///
    /// \param  query  The query sequence in which the hits were found
    /// \param  chain  A seed-set of matching positions between the query and references
    /// \return  float
    float operator()(const TString& query, const THit& hit)
    {
        float score = defaultScore_;

        // Count the number of seeds in the chain
        const size_t& referenceIdx = hit.first;
        const TSeedChain& chain = hit.second;
        size_t chainLength = seqan::length(chain);

        if (chainLength <= 0) {
            // If there are no seeds in this chain, we shouldn't be here
            return defaultScore_;
        } else if (chainLength == 1) {
            // If we have 1 seed and we can score it, return that value
            if (seedScorer_(query, chain[0], referenceIdx, score)) {
                return score;
            } else {
                return defaultScore_;
            }
        } else  // If the chain has more than 1-seed, score it differently
        {
            /*std::cout << 0 << " " 
                      << seqan::beginPositionH(chain[0]) << " "
                      << seqan::endPositionH(chain[0]) << " "
                      << seqan::endPositionH(chain[0]) - seqan::beginPositionH(chain[0]) << " "
                      //<< seqan::beginPositionH(seqan::back(chain)) << " "
                      //<< seqan::endPositionH(seqan::back(chain)) << " "
                      << CountOccurrences(query, chain[0], referenceIdx) << std::endl;*/
            // If we can score the first seed, use that as the base value
            if (seedScorer_(query, chain[0], referenceIdx, score)) {
                // Iterate over each seed in the chain past the first
                for (size_t i = 1; i < chainLength; ++i) {
                    /*std::cout << i << " " 
                              << seqan::beginPositionH(chain[i]) << " "
                              << seqan::endPositionH(chain[i]) << " "
                              << seqan::endPositionH(chain[i]) - seqan::beginPositionH(chain[i]) << " "
                              //<< seqan::beginPositionH(seqan::back(chain)) << " "
                              //<< seqan::endPositionH(seqan::back(chain)) << " "
                              << CountOccurrences(query, chain[i], referenceIdx) << std::endl;*/
                    float seedFrequency = GetFrequency(query, chain[i], referenceIdx);
                    score += std::log(seedFrequency);
                }
            } else {
                return defaultScore_;
            }
        }
        // Return the calculated pValue
        return score;
    }

    /// Calculate the expected frequency of a Kmer in
    ///
    /// Equivalent to BLASR's qLambda calcultion from LISPValueImpl.hpp
    ///
    /// \param  query  The sequence from which a substring will be counted
    /// \param  seed  A seed with the start-position of the Kmer
    /// \return  float  The frequency of the occurrence of the seed in the
    ///                 reference
    float GetFrequency(const TString& query, const TSeed& seed, const size_t& referenceIdx)
    {
        return CountOccurrences(query, seed, referenceIdx) / referenceSize_;
    }

    /// \param  query  The sequence from which a substring will be counted
    /// \param  seed  A seed with the start-position of the Kmer start position
    /// \return  the number of occurences
    size_t CountOccurrences(const TString& query, const TSeed& seed, const size_t& referenceIdx)
    {
        using namespace seqan;

        size_t occurrences = 0;
        if (LengthH(seed) == Size) {
            hash(shape_, begin(query) + beginPositionH(seed));
            for (const auto o : getOccurrences(index_, shape_))
                if (getValueI1(o) == referenceIdx) occurrences += 1;
            return occurrences;
        } else {
            hash(shape_, begin(query) + beginPositionH(seed));
            for (const auto o : getOccurrences(index_, shape_))
                if (getValueI1(o) == referenceIdx) occurrences += 1;
            return occurrences;
        }
    }

    size_t LengthH(const TSeed& seed)
    {
        using namespace seqan;

        return endPositionH(seed) - beginPositionH(seed);
    }

private:  // data
    TIndex index_;
    TShape shape_;
    size_t kmerSize_;
    float referenceSize_;
    SeedScorer<TConfig> seedScorer_;
    float defaultScore_ = 1.0f;
};
}
}  // ::PacBio::Align
