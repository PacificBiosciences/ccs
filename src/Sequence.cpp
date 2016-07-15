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

#include <stdexcept>
#include <string>

#include <pacbio/consensus/Sequence.h>

namespace PacBio {
namespace Consensus {

char Complement(char base)
{
    switch (base) {
        case 'A':
            return 'T';
        case 'a':
            return 't';
        case 'C':
            return 'G';
        case 'c':
            return 'g';
        case 'G':
            return 'C';
        case 'g':
            return 'c';
        case 'T':
            return 'A';
        case 't':
            return 'a';
        case '-':
            return '-';
        default:
            throw std::invalid_argument("invalid base!");
    }
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

}  // namespace Consensus
}  // namespace PacBio
