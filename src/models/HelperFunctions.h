// Copyright (c) 2017, Pacific Biosciences of California, Inc.
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

#include <algorithm>
#include <array>
#include <cstdint>
#include <stdexcept>
#include <utility>

#include <pacbio/consensus/ModelConfig.h>
#include <pacbio/data/Read.h>

namespace PacBio {
namespace Consensus {
namespace {

template <typename T>
inline T clip(const T val, const T (&range)[2])
{
    return std::max(range[0], std::min(val, range[1]));
}

inline uint8_t EncodeBase(const char base)
{
    const uint8_t em = detail::TranslationTable[static_cast<uint8_t>(base)];
    if (em > 3U) throw std::invalid_argument("invalid character in read!");
    return em;
}

inline uint8_t EncodeBase(const char base, const uint8_t raw_pw)
{
    if (raw_pw < 1U) throw std::runtime_error("invalid PulseWidth in read!");
    const uint8_t pw = std::min(2, raw_pw - 1);
    const uint8_t bp = detail::TranslationTable[static_cast<uint8_t>(base)];
    if (bp > 3U) throw std::invalid_argument("invalid character in read!");
    const uint8_t em = (pw << 2) | bp;
    if (em > 11U) throw std::runtime_error("read encoding error!");
    return em;
}

// first: base
// second: pw
inline std::pair<char, uint8_t> DecodeEmission(const uint8_t em)
{
    constexpr static std::array<char, 4> bases{{'A', 'C', 'G', 'T'}};
    if (em > 11U) throw std::runtime_error("encoded emission value is invalid!");
    const uint8_t base = em & 3;
    const uint8_t pw = (em >> 2) + 1;
    if (pw > 3U) throw std::runtime_error("invalid generated PulseWidth!");
    return {bases[base], pw};
}

}  // namespace anonymous
}  // namespace Consensus
}  // namespace PacBio
