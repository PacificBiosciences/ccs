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

#include <cassert>
#include <fstream>
#include <stdexcept>

#include <ssw_cpp.h>

#include <pbbam/Cigar.h>

#include <pacbio/align/SimdAlignment.h>

namespace PacBio {
namespace Align {
PariwiseAlignmentFasta SimdNeedleWunschAlignment(const std::string& target,
                                                 const std::string& query)
{
    StripedSmithWaterman::Aligner aligner;
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment;
    aligner.Align(query.c_str(), target.c_str(), target.size(), filter, &alignment);

    size_t qryPos = 0;
    size_t tgtPos = 0;

    std::string refAlign;
    std::string qryAlign;
    std::string transcript;

    for (int i = 0; i < alignment.ref_begin; ++i) {
        refAlign += target.at(tgtPos++);
        qryAlign += "-";
        transcript += "P";
    }
    for (const auto& c : BAM::Cigar::FromStdString(alignment.cigar_string)) {
        for (size_t i = 0; i < c.Length(); ++i) {
            transcript += c.Char();

            switch (c.Char()) {
                case '=':
                    assert(target.at(tgtPos) == query.at(qryPos));
                case 'M':
                case 'X':
                    refAlign += target.at(tgtPos++);
                    qryAlign += query.at(qryPos++);
                    break;
                case 'D':
                    refAlign += target.at(tgtPos++);
                    qryAlign += "-";
                    break;
                case 'I':
                case 'S':
                    refAlign += "-";
                    qryAlign += query.at(qryPos++);
                    break;
                case 'H':
                    throw std::runtime_error("H");
                default:
                    throw std::runtime_error("unknown");
            }
        }
    }

    while (tgtPos < target.size()) {
        refAlign += target.at(tgtPos++);
        qryAlign += "-";
        transcript += "P";
    }

    PariwiseAlignmentFasta result;
    result.Target = std::move(refAlign);
    result.Query = std::move(qryAlign);
    result.Transcript = std::move(transcript);

    assert(refAlign.size() == qryAlign.size());

    // std::ofstream out("out.msa");
    // out << ">ref" << std::endl
    //     << refAlign << std::endl
    //     << ">target" << std::endl
    //     << qryAlign << std::endl;
    // out.flush();
    return result;
}
}
}  // namespace PacBio::Align
