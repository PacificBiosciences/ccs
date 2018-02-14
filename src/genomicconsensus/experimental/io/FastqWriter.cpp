// Copyright (c) 2017-2018, Pacific Biosciences of California, Inc.
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

// Author: Derek Barnett

#include <pacbio/genomicconsensus/experimental/io/FastqWriter.h>

#include <stdexcept>

#include <pbbam/QualityValues.h>

using namespace std::literals::string_literals;

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

FastqWriter::FastqWriter(const Settings& settings)
    : file_{settings.fastqFilename}, out_{file_.tempFilename_}
{
    //  targetFilename, targetFilename + ".tmp"
    if (!out_) {
        throw std::runtime_error("could not open "s + file_.targetFilename_ + " for writing"s);
    }
}

void FastqWriter::Write(const std::string& header, const std::string& sequence,
                        const std::vector<uint8_t>& qualities)
{
    WriteLine('@' + header);
    WriteLine(sequence);
    WriteLine("+");
    WriteLine(PacBio::BAM::QualityValues(qualities).Fastq());
}

void FastqWriter::WriteLine(const std::string& line) { out_ << line << '\n'; }

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
