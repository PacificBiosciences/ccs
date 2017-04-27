// Copyright (c) 2016, Pacific Biosciences of California, Inc.
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

// Author: Derek Barnett

//
// Alignment blocks for banded chain alignment
//

#pragma once

#include <cfloat>

#include <vector>

#include <pbcopper/align/Seed.h>

#include <pbbam/Cigar.h>

#include <pacbio/align/BandedChainAlignment.h>

namespace PacBio {
namespace Align {
namespace Internal {

///
/// \brief The AlignInput struct acts as a "parameter-pack" for simplifiying
///        method signatures.
///
struct AlignInput
{
    const char* target;
    size_t tLen;
    const char* query;
    size_t qLen;
};

///
/// \brief The BandedGlobalAlignBlock class provides a reusable alignment
///        matrix for performing a banded, global alignment.
///
/// \note Currently only intended for use within the BandeChainAlign algorithm.
///
class BandedGlobalAlignBlock
{
public:
    BandedGlobalAlignBlock(const BandedChainAlignConfig& config) : config_(config) {}
    ~BandedGlobalAlignBlock(void) {}

public:
    ///
    /// \brief Align
    ///
    /// Aligns query to target, using a banded-global alignment, with affine
    /// gap-penalties
    ///
    /// \param target   target sequence
    /// \param tLen     target sequence length
    /// \param query    query sequence
    /// \param qLen     query sequence length
    /// \param seed     hit region
    /// \return
    ///
    PacBio::BAM::Cigar Align(const char* target, const size_t tLen, const char* query,
                             const size_t qLen, PacBio::Align::Seed seed);

private:
    std::pair<size_t, size_t> BacktraceStart(const size_t tLen, const size_t qLen) const;

    size_t IndexFor(const size_t i, const size_t j) const;

    void Init(const size_t tLen, const size_t qLen);
    size_t InitLookup(const size_t tLen, const size_t qLen);
    void InitScores(const size_t tLen, const size_t qLen, const size_t n);

private:
    const BandedChainAlignConfig& config_;

    ///
    /// \brief The LookupElement struct helps provide mappings from the 2-D
    ///        (conceptual) matrix coordinates to the linear storage vectors.
    ///
    /// \sa IndexFor
    ///
    struct LookupElement
    {
        size_t arrayStart_;
        size_t jBegin_;
        size_t jEnd_;

        LookupElement(const size_t start, const size_t jBeg, const size_t jEnd)
            : arrayStart_(start), jBegin_(jBeg), jEnd_(jEnd)
        {
        }

        LookupElement(const LookupElement&) = default;
        LookupElement(LookupElement&&) = default;
    };

    std::vector<LookupElement> lookup_;

    std::vector<float> matchScores_;
    std::vector<float> gapScores_;
};

///
/// \brief The StandardGlobalAlignBlock class probides a reusable alignment
///        matrix for standard (non-banded) global alignment.
///
/// \note Currently only intended for use within the BandeChainAlign algorithm.
///
class StandardGlobalAlignBlock
{
public:
    StandardGlobalAlignBlock(const BandedChainAlignConfig& config) : config_(config) {}
    ~StandardGlobalAlignBlock(void) {}

public:
    ///
    /// \brief Align
    ///
    /// Aligns query to target, using a global (NW) alignment, with affine
    /// gap-penalties
    ///
    /// \param target
    /// \param tLen
    /// \param query
    /// \param qLen
    /// \return
    ///
    PacBio::BAM::Cigar Align(const char* target, const size_t tLem, const char* query,
                             const size_t qLen);

private:
    std::pair<size_t, size_t> BacktraceStart(const size_t tLen, const size_t qLen) const;

    void Init(const size_t tLen, const size_t qLen);

private:
    const BandedChainAlignConfig& config_;
    typedef std::vector<std::vector<float>> Matrix_t;
    Matrix_t matchScores_;
    Matrix_t gapScores_;
};

}  // namespace Internal
}  // namespace Align
}  // namespace PacBio
