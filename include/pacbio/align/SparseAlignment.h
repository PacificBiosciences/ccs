// Author: Brett Bowman

#pragma once

#include <vector>

#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/sequence.h>

#include <pbcopper/align/Seeds.h>

#include <pacbio/align/ChainSeeds.h>
#include <pacbio/align/FindSeeds.h>
#include <pacbio/data/Sequence.h>

namespace PacBio {
namespace Align {

namespace {

///// \brief Generate an full alignment for two sequences, given their shared seeds
/////
///// \tparam    TAlignConfig  The style global/local alignment to perform
///// \tparam    TScoring      The alignment scoring scheme to use
/////
///// \param[in] seq1          The query sequence as a DnaString
///// \param[in] seq2          The reference sequence as a DnaString
///// \param[in] seeds         The seeds shared by the query and reference
///// \param[in] scoring       The alignment scoring configuration
///// \param[in] config        The alignment type configuration
/////
///// \returns  The full banded alignment between the two sequences
/////
//template <typename TAlignConfig, typename TScoring>
//seqan::Align<seqan::DnaString, seqan::ArrayGaps> SeedsToAlignment(
//    const seqan::DnaString& seq1,
//    const seqan::DnaString& seq2,
//    const seqan::SeedSet<seqan::Seed<seqan::Simple>>& seeds,
//    const TScoring& scoring,
//    const TAlignConfig& config)
//{

// TODO (dbarnett) : remove this? Doesn't look like it's used anywhere

//    using namespace seqan;

//    seqan::String<seqan::Seed<seqan::Simple>> chain;
//    chainSeedsGlobally(chain, seeds, SparseChaining());

//    seqan::Align<DnaString, ArrayGaps> alignment;
//    resize(rows(alignment), 2);
//    assignSource(row(alignment, 0), seq1);
//    assignSource(row(alignment, 1), seq2);

//    bandedChainAlignment(alignment, chain, scoring, config);

//    return alignment;
//}

/// \brief Generate an SDP alignment from two sequences
///
/// \param[in] qGramSize    qgram size to use for index hashing
/// \param[in] seq1         The query sequence
/// \param[in] seq2         The reference sequence
/// \param[in] filterHomopolymers If true, homopolymer k-mers will be filtered before searching the index.
///
/// \returns   The SDP alignment as a vector of Seeds
///
inline std::vector<Seed> SparseAlignSeeds(const size_t qGramSize, const std::string& seq1,
                                          const std::string& seq2, const bool filterHomopolymers)
{
    const auto seeds = FindSeeds(qGramSize, seq1, seq2, filterHomopolymers);
    const auto chains = ChainSeeds(seeds, ChainSeedsConfig{});
    if (chains.empty()) return std::vector<Seed>{};
    return chains[0];
}

/// \brief Generate an SDP alignment from two sequences
///
/// This overload enables homopolymer-filtering when FILTERHOMOPOLYMERS is defined.
///
/// \param[in] qGramSize    qgram size to use for index hashing
/// \param[in] seq1         The query sequence
/// \param[in] seq2         The reference sequence
///
/// \returns   The SDP alignment as a vector of Seeds
///
inline std::vector<Seed> SparseAlignSeeds(const size_t qGramSize, const std::string& seq1,
                                          const std::string& seq2)
{
#ifdef FILTERHOMOPOLYMERS
    const bool filterHomopolymers = true;
#else
    const bool filterHomopolymers = false;
#endif

    return SparseAlignSeeds(qGramSize, seq1, seq2, filterHomopolymers);
}

/// \brief Generate an SDP alignment from the best orientation of two sequences
///
/// \param[in] seq1   The query sequence as a DnaString
/// \param[in] seq2   The reference sequence as a DnaString
/// \param[in] filterHomopolymers If true, homopolymer k-mers will be filtered before searching the index.
///
/// \returns   A flag for the best orientation found, and the SDP alignment
///             from that orientation as a SeedString, in an std::pair
///
inline std::pair<size_t, std::vector<Seed>> BestSparseAlign(const std::string& seq1,
                                                            const std::string& seq2,
                                                            const bool filterHomopolymers)
{
    const auto seq2rc = ::PacBio::Data::ReverseComplement(seq2);

    const auto fwd = SparseAlignSeeds(10, seq1, seq2, filterHomopolymers);
    const auto rev = SparseAlignSeeds(10, seq1, seq2rc, filterHomopolymers);

    if (fwd.size() > rev.size()) return std::make_pair(0, fwd);
    return std::make_pair(1, rev);
}

/// \brief Generate an SDP alignment from the best orientation of two sequences
///
/// This overload enables homopolymer-filtering when FILTERHOMOPOLYMERS is defined.
///
/// \param[in] seq1   The query sequence as a DnaString
/// \param[in] seq2   The reference sequence as a DnaString
///
/// \returns   A flag for the best orientation found, and the SDP alignment
///             from that orientation as a SeedString, in an std::pair
///
inline std::pair<size_t, std::vector<Seed>> BestSparseAlign(const std::string& seq1,
                                                            const std::string& seq2)
{
#ifdef FILTERHOMOPOLYMERS
    const bool filterHomopolymers = true;
#else
    const bool filterHomopolymers = false;
#endif

    return BestSparseAlign(seq1, seq2, filterHomopolymers);
}

/// \brief Generate an SDP alignment from two sequences and hide the
///         SeqAn library dependencies
///
/// \param[in] qGramSize qgram size to use for index hashing
/// \param[in] seq1      The query sequence as an std::string
/// \param[in] seq2      The reference sequence as a std::string
/// \param[in] filterHomopolymers If true, homopolymer k-mers will be filtered before searching the index.
///
/// \returns   A vector of pairs, representing Kmer start positions
///             that match in the query and reference sequences
///
inline const std::vector<std::pair<size_t, size_t>> SparseAlign(const size_t qGramSize,
                                                                const std::string& seq1,
                                                                const std::string& seq2,
                                                                const bool filterHomopolymers)
{
    std::vector<std::pair<size_t, size_t>> result;

    const auto chain = SparseAlignSeeds(qGramSize, seq1, seq2, filterHomopolymers);
    for (const auto& s : chain)
        result.emplace_back(s.BeginPositionH(), s.BeginPositionV());

    return result;
}

/// \brief Generate an SDP alignment from two sequences and hide the
///         SeqAn library dependencies
///
/// This overload enables homopolymer-filtering when FILTERHOMOPOLYMERS is defined.
///
/// \param[in] qGramSize qgram size to use for index hashing
/// \param[in] seq1      The query sequence as an std::string
/// \param[in] seq2      The reference sequence as a std::string
///
/// \returns   A vector of pairs, representing Kmer start positions
///             that match in the query and reference sequences
///
inline const std::vector<std::pair<size_t, size_t>> SparseAlign(const size_t qGramSize,
                                                                const std::string& seq1,
                                                                const std::string& seq2)
{
#ifdef FILTERHOMOPOLYMERS
    const bool filterHomopolymers = true;
#else
    const bool filterHomopolymers = false;
#endif

    return SparseAlign(qGramSize, seq1, seq2, filterHomopolymers);
}

}  // Anonymous namespace
}
}  // ::PacBio::Align
