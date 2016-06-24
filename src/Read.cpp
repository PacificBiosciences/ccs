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

#include <cassert>

#include <pacbio/consensus/ModelConfig.h>
#include <pacbio/consensus/Read.h>

namespace PacBio {
namespace Consensus {

SNR::SNR(const double a, const double c, const double g, const double t) : A(a), C(c), G(g), T(t) {}
SNR::SNR(const std::vector<double>& snrs) : A(snrs[0]), C(snrs[1]), G(snrs[2]), T(snrs[3])
{
    assert(snrs.size() == 4);
}

namespace {
    double clamp(double val, double lo, double hi)
    {
        return std::min(std::max(val, lo), hi);
    }
}  // namespace

SNR ClampSNR(const SNR& val, const SNR& lo, const SNR& hi)
{
    return SNR(clamp(val.A, lo.A, hi.A),
               clamp(val.C, lo.C, hi.C),
               clamp(val.G, lo.G, hi.G),
               clamp(val.T, lo.T, hi.T));
}

Read::Read(const std::string& name, const std::string& seq, const std::vector<uint8_t>& ipd,
           const std::vector<uint8_t>& pw, const SNR& snr, const std::string& model)
    : Name{name}, Seq{seq}, IPD{ipd}, PulseWidth{pw}, SignalToNoise{snr}, Model{model}
{
}

MappedRead::MappedRead(const Read& read, StrandEnum strand, size_t templateStart,
                       size_t templateEnd, bool pinStart, bool pinEnd)
    : Read(read)
    , Strand{strand}
    , TemplateStart{templateStart}
    , TemplateEnd{templateEnd}
    , PinStart{pinStart}
    , PinEnd{pinEnd}
{
}

std::ostream& operator<<(std::ostream& os, const MappedRead& mr)
{
    os << "MappedRead(Read(\"" << mr.Name << "\", \"" << mr.Seq << "\", \"" << mr.Model << "\"), ";
    if (mr.Strand == StrandEnum::FORWARD)
        os << "StrandEnum_FORWARD, ";
    else if (mr.Strand == StrandEnum::REVERSE)
        os << "StrandEnum_REVERSE, ";
    else if (mr.Strand == StrandEnum::UNMAPPED)
        os << "StrandEnum_UNMAPPED, ";
    os << mr.TemplateStart << ", " << mr.TemplateEnd << ", ";
    os << mr.PinStart << ", " << mr.PinEnd << ")";
    return os;
}

}  // namespace Consensus
}  // namespace PacBio
