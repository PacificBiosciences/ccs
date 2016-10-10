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

// Author: Armin Töpfer

#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include <pbbam/BamRecord.h>

#include <pacbio/data/ArrayBase.h>

namespace PacBio {
namespace Data {
// Convert {0, 1, 2, 3, 4} to {'A', 'C', 'G', 'T', '-'}
static constexpr char TagToNucleotide(uint8_t t)
{
    switch (t) {
        case 0:
            return 'A';
        case 1:
            return 'C';
        case 2:
            return 'G';
        case 3:
            return 'T';
        case 4:
            return '-';
        default:
            return 0;
            // throw std::runtime_error("Woot is that tag? " + std::to_string(t));
    }
}
// Convert {'A', 'C', 'G', 'T', '-', 'N'} to {0, 1, 2, 3, 4, 4}
static constexpr uint8_t NucleotideToTag(char t)
{
    switch (t) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        case 'N':
            return 4;
        case '-':
            return 4;
        default:
            return 0;
            // throw std::runtime_error("Woot is that char " + std::to_string(t));
    }
}

/// A single array read that is "unrolled", as in an array of bases.
class ArrayRead
{
public:  // ctors
    ArrayRead();
    /// Constructor that needs the BamRecord to be "unrolled" and a unique index
    ArrayRead(const BAM::BamRecord& record, int idx);

    // friend std::ostream& operator<<(std::ostream& stream, const ArrayRead& r);

public:  // non-mod methods
    int ReferenceStart() const { return Record.ReferenceStart(); }
    int ReferenceEnd() const { return Record.ReferenceEnd(); }

public:  // data
    std::vector<ArrayBase> Bases;
    const BAM::BamRecord Record;
    const int Idx;
};
}  // namespace Data
}  // namespace PacBio
