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

#include <pacbio/align/SparseAlignment.h>
#include <pacbio/exception/CCSExceptions.h>

#include <pbcopper/align/Seed.h>
#include <pbcopper/qgram/Index.h>

namespace PacBio {
namespace CCS {

///
/// Find all matching seeds between two DNA sequences
///
/// \sa PacBio::Align::FindSeeds
///
/// For PacBio::CCS, homopolymer filtering is always enabled.
///
/// \param[in]  seq1    The first, or query, sequence
/// \param[in]  seq2    The second, or reference, sequence
///
/// \return Seeds collection containing all hits
///
inline PacBio::Align::Seeds FindSeeds(const size_t qGramSize, const std::string& seq1,
                                      const std::string& seq2)
{
    return PacBio::Align::FindSeeds(qGramSize, seq1, seq2, true);
}

///
/// \brief FindSeeds
///
/// \sa PacBio::Align::FindSeeds
///
/// For PacBio::CCS, homopolymer filtering is always enabled.
///
/// \param[in]  index   The hashed index on the reference sequence(s)
/// \param[in]  seq     The query sequence
/// \param[in]  qIdx    (optional) The index of the query sequence, so it can be ignored
///
/// \return map containing Seeds for each referenceIndex with a hit
///
inline std::map<size_t, PacBio::Align::Seeds> FindSeeds(const PacBio::QGram::Index& index,
                                                        const std::string& seq,
                                                        const boost::optional<size_t> qIdx)
{
    return PacBio::Align::FindSeeds(index, seq, qIdx, true);
}

///
/// \brief FindSeeds
///
/// \sa PacBio::Align::FindSeeds
///
/// For PacBio::CCS, homopolymer filtering is always enabled.
///
/// \param[in]  index   The hashed index on the reference sequence(s)
/// \param[in]  seq     The query sequence
///
/// \return map containing Seeds for each referenceIndex with a hit
///
inline std::map<size_t, PacBio::Align::Seeds> FindSeeds(const PacBio::QGram::Index& index,
                                                        const std::string& seq)
{
    return PacBio::Align::FindSeeds(index, seq, boost::none, true);
}

template <size_t TSize>
size_t CountSeeds(const PacBio::Align::Seeds& seeds)
{
    size_t count = seeds.size();
#ifdef MERGESEEDS
    for (const auto& seed : seeds)
        count += seed.Size() - TSize;
#endif
    return count;
}

template <size_t TSize>
size_t CountSeeds(const std::vector<PacBio::Align::Seed>& seeds)
{
    size_t count = seeds.size();
#ifdef MERGESEEDS
    for (const auto& seed : seeds)
        count += seed.Size() - TSize;
#endif
    return count;
}

template <size_t TSize>
size_t CountSeeds(const seqan::SeedSet<seqan::Seed<seqan::Simple>>& seeds)
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
    seqan::String<seqan::Seed<seqan::Simple>> chain;
    chainSeedsGlobally(chain, seeds, seqan::SparseChaining());

    seqan::Align<seqan::DnaString, seqan::ArrayGaps> alignment;
    seqan::resize(seqan::rows(alignment), 2);
    seqan::assignSource(seqan::row(alignment, 0), seq1);
    seqan::assignSource(seqan::row(alignment, 1), seq2);

    seqan::bandedChainAlignment(alignment, chain, scoring, config);

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

///
/// Generate an SDP alignment from two sequences
///
/// \param[in]  qGramSize   qgram size to use for index hashing
/// \param[in]  seq1        The first, or query, sequence
/// \param[in]  seq2        The second, or reference, sequence
///
/// \return The SDP alignment as a vector of Seeds
///
inline std::vector<PacBio::Align::Seed> SparseAlignSeeds(const size_t qGramSize,
                                                         const std::string& seq1,
                                                         const std::string& seq2)
{
    const auto config = PacBio::Align::ChainSeedsConfig{1, 1, 3, -1, -1, -1, INT_MAX};
    const auto seeds = PacBio::Align::FindSeeds(qGramSize, seq1, seq2, true);
    const auto chains = PacBio::Align::ChainSeeds(seeds, config);
    if (chains.empty()) return std::vector<PacBio::Align::Seed>{};
    return chains[0];
}

///
/// Generate an SDP alignment from two sequences
///
/// \param[in]  qGramSize   qgram size to use for index hashing
/// \param[in]  seq1        The first, or query, sequence
/// \param[in]  seq2        The second, or reference, sequence
/// \param seq2
///
/// \return A vector of pairs, representing Kmer start positions
///             that match in the query and reference sequences
///
inline std::vector<std::pair<size_t, size_t>> SparseAlign(const size_t qGramSize,
                                                          const std::string& seq1,
                                                          const std::string& seq2)
{
    std::vector<std::pair<size_t, size_t>> result;
    const auto chain = SparseAlignSeeds(qGramSize, seq1, seq2);
    for (const auto& s : chain)
        result.emplace_back(s.BeginPositionH(), s.BeginPositionV());
    return result;
}

}  // namespace CCS
}  // namespace PacBio
