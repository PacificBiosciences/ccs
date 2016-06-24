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

#pragma once

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace PacBio {
namespace Consensus {

struct SNR
{
    double A;
    double C;
    double G;
    double T;

    SNR(double a, double c, double g, double t);
    SNR(const std::vector<double>& snrs);

    inline double operator[](const size_t i) const
    {
        if (i == 0) return A;
        if (i == 1) return C;
        if (i == 2) return G;
        if (i == 3) return T;
        throw std::invalid_argument("SNR out of bounds!");
    }

    inline bool operator==(const SNR& other) const
    {
        return other.A == A && other.C == C && other.G == G && other.T == T;
    }

    inline bool operator!=(const SNR& other) const { return !(*this == other); }
};

SNR ClampSNR(const SNR& val, const SNR& min, const SNR& max);

struct Read
{
    Read(const std::string& name, const std::string& seq, const std::vector<uint8_t>& ipd,
         const std::vector<uint8_t>& pw, const SNR& snr, const std::string& model);
    Read(const Read& read) = default;
    Read(Read&& read) = default;

    std::string Name;
    std::string Seq;
    std::vector<uint8_t> IPD;
    std::vector<uint8_t> PulseWidth;
    SNR SignalToNoise;
    std::string Model;

    inline size_t Length() const { return Seq.length(); }
};

enum struct StrandEnum : uint8_t
{
    FORWARD,
    REVERSE,
    UNMAPPED
};

struct MappedRead : public Read
{
    MappedRead(const Read& read, StrandEnum strand, size_t templateStart, size_t templateEnd,
               bool pinStart = false, bool pinEnd = false);
    MappedRead(const MappedRead& read) = default;
    MappedRead(MappedRead&& read) = default;

    StrandEnum Strand;
    size_t TemplateStart;
    size_t TemplateEnd;
    bool PinStart;
    bool PinEnd;
};

std::ostream& operator<<(std::ostream&, const MappedRead&);

}  // namespace Consensus
}  // namespace PacBio
