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

#pragma once

#include <map>
#include <queue>
#include <utility>
#include <vector>

#include <boost/optional.hpp>

#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/sequence.h>

#include <pacbio/ccs/ChainSeeds.h>
#include <pacbio/ccs/Exceptions.h>

namespace PacBio {
namespace CCS {
namespace {

inline size_t SafeSubtract(size_t size, size_t k) { return size > k ? size - k : 0; }
template <typename TShape>
class HpHasher
{
public:
    HpHasher(TShape& shape)
    {
        using namespace seqan;
        using std::string;

        const char dna[4] = {'A', 'C', 'G', 'T'};

        for (size_t i = 0; i < 4; i++) {
            DnaString s = string(length(shape), dna[i]);
            hashes[i] = hash(shape, begin(s));
        }
    }

    inline bool operator()(const unsigned h) const
    {
        if (h == hashes[0] || h == hashes[1] || h == hashes[2] || h == hashes[3]) return true;

        return false;
    }

private:
    unsigned hashes[4];
};

}  // anonymous namespace

template <typename TConfig>
void FindSeeds(seqan::SeedSet<seqan::Seed<seqan::Simple>>& seeds, const seqan::DnaString& seq1,
               const seqan::DnaString& seq2)
{
    using namespace seqan;

    typedef Shape<Dna, typename TConfig::ShapeType> TShape;

    Index<DnaString, typename TConfig::IndexType> index(seq1);
    TShape shape = indexShape(index);
    HpHasher<TShape> isHomopolymer(shape);
    auto start = begin(seq2);
    size_t end = SafeSubtract(length(seq2) + 1, TConfig::Size);

    hashInit(shape, start);

    for (size_t j = 0; j < end; j++) {
        if (isHomopolymer(hashNext(shape, start + j))) continue;

        auto hits = getOccurrences(index, shape);

        for (const auto& i : hits) {
            Seed<Simple> seed(i, j, TConfig::Size);

#ifdef MERGESEEDS
            if (!addSeed(seeds, seed, 0, Merge()))
#endif
            {
                addSeed(seeds, seed, Single());
            }
        }
    }
}

template <typename TConfig>
void FindSeeds(
    std::map<size_t, seqan::SeedSet<seqan::Seed<seqan::Simple>>>& seeds,
    const seqan::Index<seqan::StringSet<seqan::DnaString>, typename TConfig::IndexType>& index,
    const seqan::DnaString& seq, const boost::optional<size_t> qIdx = boost::none)
{
    using namespace seqan;
    using namespace std;

    typedef Shape<Dna, typename TConfig::ShapeType> TShape;

    TShape shape = indexShape(index);
    HpHasher<TShape> isHomopolymer(shape);
    auto start = begin(seq);
    size_t end = SafeSubtract(length(seq) + 1, TConfig::Size);

    hashInit(shape, start);

    for (size_t i = 0; i < end; i++) {
        if (isHomopolymer(hashNext(shape, start + i))) continue;

        auto hits = getOccurrences(index, shape);

        for (const auto& hit : hits) {
            size_t rIdx;

            if (qIdx && (rIdx = getValueI1(hit)) == *qIdx) continue;

            size_t j = getValueI2(hit);
            Seed<Simple> seed(i, j, TConfig::Size);

#ifdef MERGESEEDS
            if (!addSeed(seeds[rIdx], seed, 0, Merge()))
#endif
            {
                addSeed(seeds[rIdx], seed, Single());
            }
        }
    }
}

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

template <size_t TSize>
void FilterSeeds(std::map<size_t, seqan::SeedSet<seqan::Seed<seqan::Simple>>>& seeds,
                 const size_t nBest)
{
    using namespace std;

    if (seeds.size() <= nBest) return;

    // keep a priority queue of the biggest hits,
    // sorted ascendingly. Bump the least value if a new one is bigger.
    priority_queue<size_t, std::vector<size_t>, std::greater<size_t>> best;

    for (const auto& kv : seeds) {
        size_t nSeeds = CountSeeds<TSize>(kv.second);

        if (best.size() < nBest) {
            best.push(nSeeds);
        } else if (nSeeds > best.top()) {
            best.pop();
            best.push(nSeeds);
        }
    }

    size_t minSize = best.top();

    for (auto it = seeds.begin(); it != seeds.end();) {
        if (CountSeeds<TSize>(it->second) < minSize) {
            it = seeds.erase(it);
        } else {
            ++it;
        }
    }
}

template <typename TAlignConfig, typename TScoring>
seqan::Align<seqan::DnaString, seqan::ArrayGaps> SeedsToAlignment(
    const seqan::DnaString& seq1, const seqan::DnaString& seq2,
    const seqan::SeedSet<seqan::Seed<seqan::Simple>>& seeds, const TScoring& scoring,
    const TAlignConfig& config)
{
    using namespace seqan;

    String<Seed<Simple>> chain;
    chainSeedsGlobally(chain, seeds, SparseChaining());

    Align<DnaString, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), seq1);
    assignSource(row(alignment, 1), seq2);

    bandedChainAlignment(alignment, chain, scoring, config);

    return alignment;
}

// TODO (lhepler) : investigate default values other than 10
template <size_t TSize = 10, typename TShape = seqan::UngappedShape<TSize>,
          typename TIndex = seqan::IndexQGram<TShape>>
struct FindSeedsConfig
{
    typedef TIndex IndexType;
    typedef TShape ShapeType;
    static const size_t Size = TSize;
};

template <size_t TSize>
seqan::String<seqan::Seed<seqan::Simple>> SparseAlign(const seqan::DnaString& seq1,
                                                      const seqan::DnaString& seq2)
{
    using namespace seqan;

    typedef FindSeedsConfig<TSize> TConfig;

    SeedSet<Seed<Simple>> seeds;
    FindSeeds<TConfig>(seeds, seq1, seq2);

    String<Seed<Simple>> chain;
    ChainSeeds(chain, seeds);

    return chain;
}

template <size_t TSize>
std::vector<std::pair<size_t, size_t>> SparseAlign(const std::string& seq1, const std::string& seq2)
{
    using namespace seqan;
    using namespace std;

    DnaString a = seq1, b = seq2;
    auto chain = SparseAlign<TSize>(a, b);

    vector<pair<size_t, size_t>> result;

    for (const auto& s : chain) {
        result.push_back(make_pair(beginPositionH(s), beginPositionV(s)));
    }

    return result;
}

}  // namespace CCS
}  // namespace PacBio
