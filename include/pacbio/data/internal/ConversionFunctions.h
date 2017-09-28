// Copyright (c) 2011-2017, Pacific Biosciences of California, Inc.
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

// Author: David Seifert

#pragma once

#include <array>
#include <cassert>
#include <cstdint>
#include <stdexcept>
#include <string>

using namespace std::literals::string_literals;  // for std::operator ""s

namespace PacBio {
namespace Data {
namespace detail {
namespace {

// IMPORTANT NOTICE
// It should be emphasized here again, that all functions in this file
// are an *implementation detail* of Unanimity. Unanimity does not
// leak ambiguous bases across public API boundaries. You cannot rely
// on any interfaces in this file when consuming Unanimity. You have
// been warned.

// 1. ASCII <-> NCBI2na
inline constexpr uint8_t ASCIIToNCBI2naImpl(const char base)
{
    // We also allow converting lowercase ASCII (a/c/g/t)
    // into their respective NCBI2na.
    constexpr const std::array<uint8_t, 256> lookupTable{
        {/*   0 -  15: */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
         /*  16 -  31: */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
         /*  32 -  47: */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
         /*  48 -  63: */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
         /*  64 -  79: */ 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
         /*  80 -  95: */ 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
         /*  96 - 111: */ 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
         /* 112 - 127: */ 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,

         /* 128 - 143: */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
         /* 144 - 159: */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
         /* 160 - 175: */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
         /* 176 - 191: */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
         /* 192 - 207: */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
         /* 208 - 223: */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
         /* 224 - 239: */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
         /* 240 - 255: */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}};

    // On a two's complement architecture, the static_cast is a no-op
    const uint8_t result = lookupTable[static_cast<uint8_t>(base)];

    // It is the user's responsibility to check whether the result
    // is 4 and then error out on that.
    return result;
}

inline constexpr char NCBI2naToASCIIImpl(const uint8_t NCBI2naBase)
{
    // 4 and higher require 3 bits, too many bits for our representation
    assert(NCBI2naBase < 4);

    constexpr const std::array<char, 4> lookupTable{{'A', 'C', 'G', 'T'}};

    return lookupTable[NCBI2naBase];
}

// 2. ASCII <-> NCBI4na
inline constexpr uint8_t ASCIIToNCBI4naImpl(const char base, const bool checkValid)
{
    constexpr const std::array<uint8_t, 256> lookupTable{
        {/*   0 -  15: */ 0, 0, 0,  0, 0,  0, 0, 0, 0,  0,  0, 0,  0, 0, 0,  0,
         /*  16 -  31: */ 0, 0, 0,  0, 0,  0, 0, 0, 0,  0,  0, 0,  0, 0, 0,  0,
         /*  32 -  47: */ 0, 0, 0,  0, 0,  0, 0, 0, 0,  0,  0, 0,  0, 0, 0,  0,
         /*  48 -  63: */ 0, 0, 0,  0, 0,  0, 0, 0, 0,  0,  0, 0,  0, 0, 0,  0,
         /*  64 -  79: */ 0, 1, 14, 2, 13, 0, 0, 4, 11, 0,  0, 12, 0, 3, 15, 0,
         /*  80 -  95: */ 0, 0, 5,  6, 8,  0, 7, 9, 0,  10, 0, 0,  0, 0, 0,  0,
         /*  96 - 111: */ 0, 1, 14, 2, 13, 0, 0, 4, 11, 0,  0, 12, 0, 3, 15, 0,
         /* 112 - 127: */ 0, 0, 5,  6, 8,  0, 7, 9, 0,  10, 0, 0,  0, 0, 0,  0,

         /* 128 - 143: */ 0, 0, 0,  0, 0,  0, 0, 0, 0,  0,  0, 0,  0, 0, 0,  0,
         /* 144 - 159: */ 0, 0, 0,  0, 0,  0, 0, 0, 0,  0,  0, 0,  0, 0, 0,  0,
         /* 160 - 175: */ 0, 0, 0,  0, 0,  0, 0, 0, 0,  0,  0, 0,  0, 0, 0,  0,
         /* 176 - 191: */ 0, 0, 0,  0, 0,  0, 0, 0, 0,  0,  0, 0,  0, 0, 0,  0,
         /* 192 - 207: */ 0, 0, 0,  0, 0,  0, 0, 0, 0,  0,  0, 0,  0, 0, 0,  0,
         /* 208 - 223: */ 0, 0, 0,  0, 0,  0, 0, 0, 0,  0,  0, 0,  0, 0, 0,  0,
         /* 224 - 239: */ 0, 0, 0,  0, 0,  0, 0, 0, 0,  0,  0, 0,  0, 0, 0,  0,
         /* 240 - 255: */ 0, 0, 0,  0, 0,  0, 0, 0, 0,  0,  0, 0,  0, 0, 0,  0}};

    const uint8_t result = lookupTable[static_cast<uint8_t>(base)];

    if ((result == 0) && (checkValid))
        throw std::runtime_error("Invalid ASCII value ('"s + base + "', ordinal "s +
                                 std::to_string(static_cast<int>(base)) +
                                 ") for converting into NCBI4na format!"s);

    return result;
}

inline constexpr char NCBI4naToASCIIImpl(const uint8_t NCBI4naBase)
{
    // NCBI4na 0, i.e., gaps are non-sensical for our use-cases
    assert(NCBI4naBase != 0);

    // 16 and higher require 5 bits, too many bits for our representation
    assert(NCBI4naBase < 16);

    constexpr const std::array<char, 16> lookupTable{
        {'\0', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'}};

    return lookupTable[NCBI4naBase];
}

// 3. NCBI2na <-> NCBI4na
inline constexpr uint8_t NCBI2naToNCBI4naImpl(const uint8_t NCBI2naBase)
{
    assert(NCBI2naBase < 4);

    const uint8_t result = (1 << NCBI2naBase);
    return result;
}

inline constexpr uint8_t NCBI4naToNCBI2naImpl(const uint8_t NCBI4naBase)
{
    assert(NCBI4naBase < 16);

    // 4 represents an invalid sentinel value
    constexpr const std::array<uint8_t, 16> lookupTable{{
        /* -/0  */ 4, /* A/1  */ 0, /* C/2  */ 1, /* M/3  */ 4,
        /* G/4  */ 2, /* R/5  */ 4, /* S/6  */ 4, /* V/7  */ 4,
        /* T/8  */ 3, /* W/9  */ 4, /* Y/10 */ 4, /* H/11 */ 4,
        /* K/12 */ 4, /* D/13 */ 4, /* B/14 */ 4, /* N/15 */ 4,
    }};

    const uint8_t result = lookupTable[NCBI4naBase];

    if (result > 3)
        throw std::runtime_error("Invalid NCBI4na value for converting into NCBI2na format!");

    return result;
}

// 4. Utility functions
//
//    - for checking number of bits in NCBI4na,
//      i.e., the ploidy of the position
//    - for compressing a major and minor allele
//      into a single (ambiguous) base. This is a
//      stopgap solution and will be removed in
//      mid-term future.

inline constexpr uint8_t numSetBitsImpl(const uint8_t NCBI4naBase)
{
    assert(NCBI4naBase < 16);

    // Also, uses lookup table again
    // http://gurmeet.net/puzzles/fast-bit-counting-routines/
    constexpr const std::array<uint8_t, 16> lookupTable{
        {0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4}};
    return lookupTable[NCBI4naBase];
}

inline constexpr char createAmbiguousBase(const char firstBase, const char secondBase)
{
    const uint8_t firstNCBI4na = ASCIIToNCBI4naImpl(firstBase, true);
    const uint8_t secondNCBI4na = ASCIIToNCBI4naImpl(secondBase, true);

    const uint8_t combinedBaseNCBI4na = (firstNCBI4na | secondNCBI4na);
    const char result = NCBI4naToASCIIImpl(combinedBaseNCBI4na);

    return result;
}

}  // namespace
}  // namespace detail
}  // namespace Data
}  // namespace PacBio
