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
#include <map>
#include <vector>

#include <pacbio/data/ArrayRead.h>
#include <pacbio/data/FisherResult.h>

namespace PacBio {
namespace Data {
class MSAColumn
{
public:
    // Relative per nucleotide abundance
    double Frequency(int i) const { return (*this)[i] / static_cast<double>(Coverage()); }
    double Frequency(char c) const { return (*this)[c] / static_cast<double>(Coverage()); }

    // Nucleotide counts
    int operator[](int i) const { return counts[i]; }
    int& operator[](int i) { return counts[i]; }
    int operator[](char c) const { return counts[NucleotideToTag(c)]; }
    int& operator[](char c) { return counts[NucleotideToTag(c)]; }

    operator std::array<int, 5>() { return counts; }
    explicit operator int() { return Coverage(); }

public:
    int Coverage() const;

public:
    void AddFisherResult(const FisherResult& f);
    void AddFisherResult(const std::map<std::string, double>& f);

public:
    std::array<int, 5> counts{{0, 0, 0, 0, 0}};
    std::map<std::string, int> insertions;
    std::map<std::string, double> insertionsPValues;
    std::array<double, 5> pValues{{1, 1, 1, 1, 1}};
    std::array<double, 5> mask{{0, 0, 0, 0, 0}};
    bool hit = false;
    int argMax = 0;

public:
    friend std::ostream& operator<<(std::ostream& stream, const MSAColumn& r)
    {
        for (int j = 0; j < 5; ++j)
            stream << r.counts.at(j) << "\t" << r.pValues.at(j) << "\t";
        return stream;
    }
};
}  // namespace Data
}  // namespace PacBio
