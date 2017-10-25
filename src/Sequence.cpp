// Copyright (c) 2011-2016, Pacific Biosciences of California, Inc.
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

#include <array>
#include <stdexcept>
#include <string>

using namespace std::literals::string_literals;  // for std::operator ""s

#include <pacbio/data/Sequence.h>

namespace PacBio {
namespace Data {

char Complement(const char base)
{
    constexpr const std::array<char, 256> lookupTable{
        {/*   0 -   7: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /*   8 -  15: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /*  16 -  23: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /*  24 -  31: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /*  32 -  39: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /*  40 -  47: */ 0,   0,   0,   0,   0,   '-', 0,   0,
         /*  48 -  55: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /*  56 -  63: */ 0,   0,   0,   0,   0,   0,   0,   0,

         /*  64 -  71: */ 0,   'T', 'V', 'G', 'H', 0,   0,   'C',
         /*  72 -  79: */ 'D', 0,   0,   'M', 0,   'K', 'N', 0,
         /*  80 -  87: */ 0,   0,   'Y', 'S', 'A', 0,   'B', 'W',
         /*  88 -  95: */ 0,   'R', 0,   0,   0,   0,   0,   0,

         /*  96 - 103: */ 0,   't', 'v', 'g', 'h', 0,   0,   'c',
         /* 104 - 111: */ 'd', 0,   0,   'm', 0,   'k', 'n', 0,
         /* 112 - 119: */ 0,   0,   'y', 's', 'a', 0,   'b', 'w',
         /* 120 - 127: */ 0,   'r', 0,   0,   0,   0,   0,   0,

         /* 128 - 135: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 136 - 143: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 144 - 151: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 152 - 159: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 160 - 167: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 168 - 175: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 176 - 183: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 184 - 191: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 192 - 199: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 200 - 207: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 208 - 215: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 216 - 223: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 224 - 231: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 232 - 239: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 240 - 247: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 248 - 255: */ 0,   0,   0,   0,   0,   0,   0,   0}};

    const char result = lookupTable[static_cast<uint8_t>(base)];

    if (result == 0) throw std::invalid_argument(base + " is an invalid base!"s);

    return result;
}

std::string Complement(const std::string& input)
{
    std::string output;
    output.reserve(input.length());
    for (const char b : input)
        output.push_back(Complement(b));
    return output;
}

std::string Reverse(const std::string& input)
{
    std::string output;
    output.reserve(input.length());
    for (auto it = input.crbegin(); it != input.crend(); ++it)
        output.push_back(*it);
    return output;
}

std::string ReverseComplement(const std::string& input)
{
    std::string output;
    output.reserve(input.length());
    for (auto it = input.crbegin(); it != input.crend(); ++it)
        output.push_back(Complement(*it));
    return output;
}

}  // namespace Data
}  // namespace PacBio
