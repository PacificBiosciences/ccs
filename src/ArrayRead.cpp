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

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include <pbbam/BamRecord.h>

#include <pacbio/data/ArrayBase.h>

#include <pacbio/data/ArrayRead.h>

namespace PacBio {
namespace Data {

ArrayRead::ArrayRead() : Idx(-1){};
ArrayRead::ArrayRead(const BAM::BamRecord& record, int idx)
    : Record(record), Idx(idx)  // Record(std::forward<BAM::BamRecord>(record))
{
    const auto seq = Record.Sequence(BAM::Orientation::GENOMIC, true, true);
    const auto qual = Record.Qualities(BAM::Orientation::GENOMIC, true, true);
    std::string cigar;
    cigar.reserve(seq.size());
    for (const auto c : Record.CigarData(true))
        for (size_t i = 0; i < c.Length(); ++i)
            cigar += c.Char();

    BAM::QualityValues subQV;
    BAM::QualityValues delQV;
    BAM::QualityValues insQV;

    bool richQVs = Record.HasDeletionQV() && Record.HasDeletionQV() && Record.HasInsertionQV();
    if (richQVs) {
        subQV = Record.SubstitutionQV(BAM::Orientation::GENOMIC, true, true);
        delQV = Record.DeletionQV(BAM::Orientation::GENOMIC, true, true);
        insQV = Record.InsertionQV(BAM::Orientation::GENOMIC, true, true);
    }
    assert(cigar.size() == seq.size() && seq.size() == qual.size());

    Bases.reserve(cigar.length());
    if (richQVs)
        for (size_t i = 0; i < cigar.length(); ++i)
            Bases.emplace_back(cigar.at(i), seq.at(i), qual.at(i), subQV.at(i), delQV.at(i),
                               insQV.at(i));
    else
        for (size_t i = 0; i < cigar.length(); ++i)
            Bases.emplace_back(cigar.at(i), seq.at(i), qual.at(i));
}

#if __cplusplus < 201402L  // C++11
char TagToNucleotide(uint8_t t)
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

uint8_t NucleotideToTag(char t)
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
#endif
}  // namespace Data
}  // namespace PacBio
