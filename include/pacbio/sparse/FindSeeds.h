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

#include <map>
#include <vector>

#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/sequence.h>

#include <pacbio/sparse/FindSeedsConfig.h>
#include <pacbio/sparse/ChainSeeds.h>
#include <pacbio/sparse/HomopolymerHasher.h>

/* 
 * This file contains a few minimal wrapper functions around the index types
 * provided by SeqAn for finding K-mer seeds between some query sequence
 * and some reference sequence or index.
 *
 * NOTE: Though these functions should in theory work with any index type
 * supported by SeqAn, they have only been extensively tested with the QGram
 * index specialization.  Use other types at your own risk.
 * 
 * In addition to the templating around the TConfig, of which more details
 * are in FindSeedConfig.h, there are two pre-processor directives that
 * can be used to further customize the code, described below.  Speed costs
 * associated with either are highly application-specific, so we recommend
 * testing both for each new application.
 *
 * MERGESEEDS: There are two common specializations of the addSeeds function
 *      that we use to add a K-mer hit to a SeedSet - Single() and Merge().  
 *      By default FindSeeds uses Single(), which is the fastest method 
 *      because it does no chaining at all.  Alternatively, we can use the 
 *      Merge() function from SeqAn, which combines seeds that precisely 
 *      overlap with each other but is slightly slower.  The advantage of
 *      Merge() is that the resulting SeedSets are smaller, so down-stream
 *      processes that require for sorting and iterating may be greatly
 *      expedited.  Enable this directive to use Merge() instead of Single().
 *
 * FILTERHOMOPOLYMERS: By default, FindSeeds treats all K-mer seeds it finds
 *      as equal.  However if the sequences might contain large homopolymers, 
 *      or if there is a large number of sequences in the the reference / 
 *      index, it may be faster in the long run to spend some CPU cycles 
 *      checking whether a K-mer is a homopolymer before searching the index 
 *      for it.  Enable this directive to enable that filter.
 */

namespace PacBio {
namespace SparseAlignment {

/// Unsigned-integer safe subtraction - returns either the difference
/// between the arguments or zero, whichever is larger.
///
/// \param  size  The first integer, or minuend
/// \param  k  The second interger, or the subtrahend
///
/// \return  size_t  The difference
inline
size_t SafeSubtract(size_t size, size_t k)
{
    return size > k ? size - k : 0;
}

/// Find all matching seeds between two DNA sequences
///
/// \param  seeds  The SeedSet object to store the results in
/// \param  seq1  The first, or query, sequence
/// \param  seq2  The second, or reference, sequence
template<typename TConfig>
void FindSeeds(seqan::SeedSet<seqan::Seed<seqan::Simple>>* seeds,
               const seqan::DnaString& seq1, 
               const seqan::DnaString& seq2)
{
    using namespace seqan;

    typedef Shape<Dna, typename TConfig::ShapeType> TShape;

    Index<DnaString, typename TConfig::IndexType> index(seq1);
    TShape shape = indexShape(index);

#ifdef FILTERHOMOPOLYMERS
    HomopolymerHasher<TShape> isHomopolymer(shape);
#endif

    auto start = begin(seq2);
    size_t end = SafeSubtract(length(seq2) + 1, TConfig::Size);

    hashInit(shape, start);
    for (size_t i = 0; i < end; i++)
    {
        hashNext(shape, start + i);
#ifdef FILTERHOMOPOLYMERS
        if (isHomopolymer(hashNext(shape, start + i)))
            continue;
#endif

        auto hits = getOccurrences(index, shape);

        for (const auto& hit : hits)
        {
            Seed<Simple> seed(hit, i, TConfig::Size);

#ifdef MERGESEEDS
            if (!addSeed(*seeds, seed, 0, Merge()))
#endif
            {
                addSeed(*seeds, seed, Single());
            }
        }
    }
}

/// Find all matching seeds between a DNA index and the sequences 
/// represented in some supplied index of the type specified in TConfig.
/// Since some index types, most notably the QGram index, can store seeds
/// from multiple references, the return value has to be a map of seed sets 
/// rather than a single one.
///
/// \param  seeds  A map of integers-SeedSet pairs for storing results
/// \param  index  The index of of the various
/// \param  seq  The query sequence
template<typename TConfig>
void FindSeeds(std::map<size_t, seqan::SeedSet<seqan::Seed<seqan::Simple>>>* seeds,
               const seqan::Index<seqan::StringSet<seqan::DnaString>, typename TConfig::IndexType>& index,
               const seqan::DnaString& seq)
{
    using namespace seqan;
    using namespace std;

    typedef Shape<Dna, typename TConfig::ShapeType> TShape;

    TShape shape = indexShape(index);

#ifdef FILTERHOMOPOLYMERS
    HomopolymerHasher<TShape> isHomopolymer(shape);
#endif

    auto start = begin(seq);
    size_t end = SafeSubtract(length(seq) + 1, TConfig::Size);

    hashInit(shape, start);

    for (size_t i = 0; i < end; i++)
    {
        hashNext(shape, start + i);

#ifdef FILTERHOMOPOLYMERS
        if (isHomopolymer(hashNext(shape, start + i)))
            continue;
#endif

        auto hits = getOccurrences(index, shape);

        for (const auto& hit : hits)
        {
            size_t rIdx = getValueI1(hit);

            size_t j = getValueI2(hit);
            Seed<Simple> seed(i, j, TConfig::Size);

#ifdef MERGESEEDS
            if (!addSeed(seeds->operator[](rIdx), seed, 0, Merge()))
#endif
            {
                addSeed(seeds->operator[](rIdx), seed, Single());
            }
        }
    }
}

/// Find all matching seeds between a DNA index and the sequences 
/// represented in some supplied index of the type specified in TConfig.
/// Since some index types, most notably the QGram index, can store seeds
/// from multiple references, the return value has to be a map of seed sets 
/// rather than a single one.  In addition the query sequence may itself
/// be in the index, in which case we pass in it's known index so we do
/// not count it.
///
/// \param  seeds  A map of integers-SeedSet pairs for storing results
/// \param  index  The index of of the various
/// \param  seq  The query sequence
/// \param  qIdx  The index of the query in the ... index, so it can be ignored
template<typename TConfig>
void FindSeeds(std::map<size_t, seqan::SeedSet<seqan::Seed<seqan::Simple>>>* seeds,
               const seqan::Index<seqan::StringSet<seqan::DnaString>, typename TConfig::IndexType>& index,
               const seqan::DnaString& seq,
               const size_t qIdx)
{
    using namespace seqan;
    using namespace std;
    
    typedef Shape<Dna, typename TConfig::ShapeType> TShape;

    TShape shape = indexShape(index);

#ifdef FILTERHOMOPOLYMERS
    HomopolymerHasher<TShape> isHomopolymer(shape);
#endif

    auto start = begin(seq);
    size_t end = SafeSubtract(length(seq) + 1, TConfig::Size);

    hashInit(shape, start);

    for (size_t i = 0; i < end; i++)
    {
        hashNext(shape, start + i);

#ifdef FILTERHOMOPOLYMERS
        if (isHomopolymer(hashNext(shape, start + i)))
            continue;
#endif

        auto hits = getOccurrences(index, shape);

        for (const auto& hit : hits)
        {
            size_t rIdx;

            if ((rIdx = getValueI1(hit)) == qIdx)
                continue;

            size_t j = getValueI2(hit);
            Seed<Simple> seed(i, j, TConfig::Size);

#ifdef MERGESEEDS
            if (!addSeed(seeds->operator[](rIdx), seed, 0, Merge()))
#endif
            {
                addSeed(seeds->operator[](rIdx), seed, Single());
            }
        }
    }
}

} // SparseAlignment
} // PacBio
