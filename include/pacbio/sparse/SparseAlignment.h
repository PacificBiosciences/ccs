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

// Author: Brett Bowman

#pragma once

#include <vector>

#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/sequence.h>

#include <pacbio/sparse/FindSeeds.h>
#include <pacbio/sparse/FindSeedsConfig.h>
#include <pacbio/sparse/ChainSeeds.h>

namespace PacBio {
namespace SparseAlignment {

namespace {

/// \brief Generate an full alignment for two sequences, given their shared seeds
///
/// \tparam    TAlignConfig  The style global/local alignment to perform
/// \tparam    TScoring      The alignment scoring scheme to use
///
/// \param[in] seq1          The query sequence as a DnaString
/// \param[in] seq2          The reference sequence as a DnaString
/// \param[in] seeds         The seeds shared by the query and reference
/// \param[in] scoring       The alignment scoring configuration
/// \param[in] config        The alignment type configuration
///
/// \returns  The full banded alignment between the two sequences 
///
template<typename TAlignConfig, typename TScoring>
const seqan::Align<seqan::DnaString, seqan::ArrayGaps>
SeedsToAlignment(const seqan::DnaString& seq1, const seqan::DnaString& seq2,
                 const seqan::SeedSet<seqan::Seed<seqan::Simple>>& seeds,
                 const TScoring& scoring,
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

/// \brief Generate an SDP alignment from two sequences
///
/// \tparam    TSize   The size of the seeds to look for
///
/// \param[in] seq1   The query sequence as a DnaString
/// \param[in] seq2   The reference sequence as a DnaString
///
/// \returns   The SDP alignment as a SeedString
/// 
template <size_t TSize>
const seqan::String<seqan::Seed<seqan::Simple>>
SparseAlign(const seqan::DnaString& seq1, const seqan::DnaString& seq2)
{
    using namespace seqan;

    typedef FindSeedsConfig<TSize> TConfig;
    ChainSeedsConfig chainConfig;

    SeedSet<Seed<Simple>> seeds;
    FindSeeds<TConfig>(&seeds, seq1, seq2);

    std::vector<String<Seed<Simple>>> chains;
    String<Seed<Simple>> chain;
    ChainSeeds(&chains, seeds, chainConfig);

    if (chains.size() == 0)
        return chain;

    return chains[0];
}

/// \brief Generate an SDP alignment from the best orientation of two sequences
///
/// \param[in] seq1   The query sequence as a DnaString
/// \param[in] seq2   The reference sequence as a DnaString
///
/// \returns   A flag for the best orientation found, and the SDP alignment 
///             from that orientation as a SeedString, in an std::pair
/// 
const std::pair<size_t, seqan::String<seqan::Seed<seqan::Simple>>>
BestSparseAlign(const seqan::DnaString& seq1, const seqan::DnaString& seq2)
{
    using namespace seqan;

    //TODO(bbowman) Clean this up, pronto
    DnaString seq2_str(seq2); 
    reverseComplement(seq2_str);
    const std::string seq2rc(toCString(CharString(seq2_str)));

    const auto& fwd = SparseAlign<10>(seq1, seq2);
    const auto& rev = SparseAlign<10>(seq1, seq2rc);

    if (length(fwd) > length(rev))
        return std::make_pair(0, fwd);
    return std::make_pair(1, rev);
}

/// \brief Generate an SDP alignment from two sequences and hide the 
///         SeqAn library dependencies
///
/// \param[in] seq1   The query sequence as an std::string
/// \param[in] seq2   The reference sequence as a std::string
///
/// \returns   A vector of pairs, representing Kmer start positions
///             that match in the query and reference sequences
/// 
template <size_t TSize>
const std::vector<std::pair<size_t, size_t>>
SparseAlign(const std::string& seq1, const std::string& seq2)
{
    using namespace seqan;
    using namespace std;

    DnaString a = seq1, b = seq2;
    auto chain = SparseAlign<TSize>(a, b);

    vector<pair<size_t, size_t>> result;

    for (const auto& s : chain)
    {
        result.push_back(make_pair(beginPositionH(s), beginPositionV(s)));
    }

    return result;
}

}  // Anonymous namespace
}}  // ::PacBio::SparseAlignment
