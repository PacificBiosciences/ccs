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

// Programmer: Brett Bowman

#pragma once

#include <stdbool.h>

#include <seqan/seeds.h>
#include <seqan/index.h>

namespace SparseAlignment {

template<typename TConfig>
class SeedScorer
{

public:  // Template types and values
   typedef typename TConfig::IndexType IndexType;
   typedef typename TConfig::ShapeType ShapeType;
   size_t Size = TConfig::Size;
   
public:  // Derived types and values
   typedef seqan::Index<seqan::StringSet<seqan::DnaString>, IndexType> TIndex;
   typedef seqan::Finder<TIndex> TFinder;
   typedef seqan::Dna TAlphabet;
   typedef typename seqan::ValueSize<TAlphabet>::Type TAlphabetSize;
   // Equivalent to BLASR's "T_Tuple"
   typedef seqan::Shape<seqan::Dna, ShapeType> TShape;
   // Equivalent to BLASR's "MatchPos"
   typedef seqan::Seed<seqan::Simple> TSeed;

public:  // structors
    SeedScorer(const TIndex& index, 
               const size_t kmerSize)
        : index_(index)
        , shape_(seqan::indexShape(index))
        , kmerFinder_(index)
        , kmerSize_(kmerSize)
        , alphabetSize_(seqan::ValueSize<seqan::Dna>::VALUE)
    {}

    // Move constructor
    SeedScorer(SeedScorer&& src) = delete;
    // Copy constructor
    SeedScorer(const SeedScorer& src) = delete;
    // Move assignment operator
    SeedScorer& operator=(SeedScorer&& rhs) = delete;
    // Copy assignment operator
    SeedScorer& operator=(const SeedScorer& rhs) = delete;
    // Destructor
    ~SeedScorer() = default;

public:  // non-modifying methods

    /// Score a given seed found that matches the reference index, returning
    /// a number approximating it's log-likehood.
    ///
    /// Equivalent to BLASR's PMatch function.
    ///
    /// \param  seed  The position and length of the query
    /// \param  query  The query sequence being
    /// \param  pValue  The number of occurences
    /// \return  bool  Flag whether the scoring succeeded
    bool operator()(const seqan::DnaString& query,
                    const TSeed& seed,
                    const size_t referenceIdx,
                    float& score)
    {
        // Default
        score = 1.0f;

        // If the count is 0, something when wrong - return failure;
        size_t count = CountOccurrences(query, seed, referenceIdx);
        //std::cout << beginPositionH(seed) << " " << endPositionH(seed) << " " << count << std::endl;
        //std::cout << seqan::infix(query, beginPositionH(seed), endPositionH(seed)) << std::endl;
        if (count == 0) return false;

        // If the match is shorter than our Kmer-size, just assume we will
        // always find a match (e.g. return 0, since it's Log(1.0))
        size_t currSeedSize = seqan::seedSize(seed);
        if (currSeedSize < kmerSize_)
        {
            score = 0.0f;
            return true;
        }
       
        /*std::cout << "S " << 0 << " " 
                  << seqan::beginPositionV(seed) << " "
                  << seqan::endPositionV(seed) << " "
                  << seqan::endPositionV(seed) - seqan::beginPositionH(seed) << " "
                  << count << " "
                  << score << std::endl;*/
                  
        // Iterate over each successor seed in the match
        for (size_t i = 1; i <= currSeedSize - kmerSize_; ++i)
        {
            //std::cout << "i:" << i << " " << seqan::beginPositionH(seed) << " " << currSeedSize << " " << Size << std::endl;
            size_t nextShapeCount = 0;
            size_t nextPossibleCount = CountPossibleSuccessors(query, seed, i, referenceIdx, nextShapeCount);

            /*std::cout << "S " << i << " " 
                      << seqan::beginPositionV(seed) << " "
                      << seqan::endPositionV(seed) << " "
                      << seqan::endPositionV(seed) - seqan::beginPositionV(seed) << " "
                      //<< CountOccurrences(query, seed, referenceIdx) << " "
                      << nextShapeCount << " " 
                      << nextPossibleCount << " "
                      << score << std::endl;*/

            // There is no background distribution available for this
            // sequence context, therefore no way to evaluate p-value.
            if (nextPossibleCount == 0) 
                return false;

            // Add the log-probability of seeing this particule successor to
            // the running sum of the pValue
            score += std::log((nextShapeCount / static_cast<float>(nextPossibleCount)));
        }

        // Return success
        return true;
    }

    /// Count the number of times the possible successor Kmers to a given
    /// seed are found,
    ///
    /// Equivalent to BLASR's SumRightShiftMarginalTupleCounts function.
    ///
    /// \param  TSeed  the position and length of the query
    /// \return  size_t  the total occurrences of all possible successors
    size_t CountPossibleSuccessors(const seqan::DnaString& query,
                                   const TSeed& seed,
                                   const size_t offset,
                                   const size_t referenceIdx,
                                   size_t& nextShapeCount)
    {
        // Initialize return values
        nextShapeCount = 0;
        size_t nextPossibleCount = 0;

        // Calculate the relative
        size_t beginPos = seqan::beginPositionH(seed) + offset;
        size_t endPos = beginPos + kmerSize_;

        // Pull-out pointers to the successor string and it's N-1 root
        seqan::DnaString referenceString = seqan::infix(query, beginPos, endPos);
        seqan::DnaString rootString = seqan::prefix(referenceString, kmerSize_-1);
        //std::cout << referenceString << " " 
        //          << rootString << std::endl;

        // Count the number of occurrences of all possible successor strings
        for (TAlphabetSize i = 0; i < alphabetSize_; ++i)
        {
            seqan::DnaString testString = rootString;
            seqan::append(testString, alphabet_[i]);
            size_t testStringCount = CountOccurrences(testString, referenceIdx);
            //std::cout << referenceString << " " 
            //          << rootString << " " 
            //          << alphabet_[i] << " " 
            //          << testString << " "
            //          << testStringCount << std::endl;

            // Add the count for the current testString to the running total
            nextPossibleCount += testStringCount;

            // If the test successor is equal to the actual successor, save it.
            if (testString == referenceString)
                nextShapeCount = testStringCount;
        }

        // Return the sum of the counts of all possible successor strings
        // the count for the actual successor is returned via reference.
        return nextPossibleCount;
    }

    /// Count the number of times a given sequence is found in the index, 
    /// usually representing all Kmers found in a reference.
    ///
    /// Equivalent to BLASR's GetTupleCount function.
    ///
    /// \param  query  The sequence whose occurrences are to be counted
    /// \return  the number of occurences
    size_t CountOccurrences(const seqan::DnaString& query,
                            const size_t referenceIdx)
    {
        using namespace seqan;

        //std::cout << "CO " << query << " " << referenceIdx << std::endl;
        size_t occurrences = 0;
        clear(kmerFinder_);
        while (find(kmerFinder_, query))
        {
            //std::cout << position(kmerFinder_) << std::endl;
            if (getValueI1(position(kmerFinder_)) == referenceIdx)
                occurrences += 1;
        }
        //std::cout << "CO " << query << " " << referenceIdx << " " << occurrences << std::endl;
        return occurrences;
    }

    /// Count the number of times a given subtring sequence is found in the
    /// index, usually representing all Kmers found in a reference.
    ///
    /// Equivalent to BLASR's GetTupleCount function.
    ///
    /// \param  query  The sequence from which a substring will be counted
    /// \param  seed  A seed with the start-position of the Kmer start position
    /// \return  the number of occurences
    size_t CountOccurrences(const seqan::DnaString& query,
                            const TSeed& seed,
                            const size_t referenceIdx)
    {
        using namespace seqan;

        size_t occurrences = 0;
        hash(shape_, begin(query) + beginPositionH(seed));
        for (const auto o : getOccurrences(index_, shape_))
        {
            if (getValueI1(o) == referenceIdx)
                occurrences += 1;
        }
        return occurrences;
    }

private:  // data
    TIndex index_;
    TShape shape_;
    TFinder kmerFinder_;
    size_t kmerSize_;
    TAlphabetSize alphabetSize_;
    char alphabet_ [4] = {'A', 'C', 'G', 'T'};

};

}  // SparseAlignment
