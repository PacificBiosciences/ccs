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

// Author: Armin Töpfer

#pragma once

#include <array>
#include <cmath>
#include <exception>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <list>
#include <memory>
#include <numeric>
#include <sstream>
#include <unordered_map>
#include <vector>

#include <pacbio/data/MSA.h>
#include <pbcopper/json/JSON.h>

#undef major  // gcc hack
namespace PacBio {
namespace Juliet {

struct VariantNucleotide
{
    VariantNucleotide(char nucleotide) : frequency(1), pValue(0), nucleotide(nucleotide), major(1)
    {
    }
    VariantNucleotide(char nucleotide, double frequency, double pValue)
        : frequency(frequency), pValue(pValue), nucleotide(nucleotide), major(0)
    {
    }
    double frequency;
    double pValue;
    char nucleotide;
    bool major;
};

/// Given a MSA and p-values for each nucleotide of each position,
/// generate machine-interpretable and human-readable output about mutated
/// amino acids.
class ResistanceCaller
{
public:
    /// Constructor needs a multiple sequence alignment.
    ResistanceCaller(const Data::MSA& msa);

public:
    /// Generate JSON output of variant amino acids
    JSON::Json JSON();

public:
    /// Generate pretty print output of variant amino acids
    static void Print(std::ostream& out, const JSON::Json& j, bool onlyKnownDRMs, bool details);

    /// Generate HTML output of variant amino acids
    static void HTML(std::ostream& out, const JSON::Json& j, bool onlyKnownDRMs, bool details);

private:
    double MaxFrequency(std::vector<VariantNucleotide> codon);

    void AddPosition(std::vector<VariantNucleotide>&& nucs);

    inline char Ref(int i) const { return ref_[i]; }

    std::string CodonRef(int hxb2Position) const;

    char AminoacidRef(int hxb2Position) const;

    inline std::string CodonString(const std::vector<VariantNucleotide>& codon) const;

    std::vector<std::vector<VariantNucleotide>> CreateCodons(const int hxb2Position) const;

private:
    Data::MSA msa_;
    std::vector<std::vector<VariantNucleotide>> nucleotides_;
    int begin_;
    int end_;

    static const std::unordered_map<int, std::string> resistantCodon_;
    static const std::unordered_map<std::string, char> codonToAmino_;
    static const std::string ref_;
};
}
}  // ::PacBio::Juliet