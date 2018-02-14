// Copyright (c) 2018, Pacific Biosciences of California, Inc.
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

#include <pacbio/genomicconsensus/experimental/io/VcfWriter.h>

#include <chrono>
#include <sstream>
#include <stdexcept>

#include <pbbam/DataSet.h>

#include <pacbio/UnanimityVersion.h>

using namespace std::literals::string_literals;

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

VcfWriter::VcfWriter(const Settings& settings, const std::vector<ReferenceWindow>& refWindows)
    : file_{settings.vcfFilename}, out_{file_.tempFilename_}
{
    //  targetFilename, targetFilename + ".tmp"
    if (!out_) {
        throw std::runtime_error("could not open "s + file_.targetFilename_ + " for writing"s);
    }

    WriteLine("##fileformat=VCFv4.3");
    WriteLine("##fileDate=" + PacBio::BAM::ToIso8601(std::chrono::system_clock::now()));
    WriteLine("##source=" + PacBio::UnanimityVersion());
    WriteLine("##reference=file://" + settings.referenceFilename);

    for (const auto& ref : refWindows)
        WriteLine("##contig=<ID=" + ref.name + ",length=" + std::to_string(ref.Length()) + ">");

    WriteLine("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
}

void VcfWriter::WriteLine(const std::string& line) { out_ << line << '\n'; }

void VcfWriter::WriteVariant(const Variant& v)
{
    // "VCF lite"
    auto pos = v.refStart;
    std::string ref;
    std::string alt;

    // indel
    if (v.refSeq.empty() || v.readSeq1.empty() || (v.IsHeterozygous() && v.readSeq2->empty())) {
        ref = v.refPrev + v.refSeq;
        alt = v.readPrev + v.readSeq1;
        if (v.IsHeterozygous()) {
            alt.push_back(',');
            alt += v.readPrev;
            alt += v.readSeq2.get();
        }
    }

    // substitution
    else {
        ++pos;
        ref = v.refSeq;
        if (v.IsHeterozygous()) {
            if (v.refSeq == v.readSeq1)
                alt = v.readSeq2.get();
            else if (v.refSeq == v.readSeq2.get())
                alt = v.readSeq1;
            else {
                alt = v.readSeq1;
                alt.push_back(',');
                alt += v.readSeq2.get();
            }
        } else
            alt = v.readSeq1;
    }

    std::stringstream line;
    line << v.refName << '\t' << pos << "\t.\t" << ref << '\t' << alt << '\t' << v.confidence.get()
         << "\tPASS";
    WriteLine(line.str());
}

void VcfWriter::WriteVariants(const std::vector<Variant>& variants)
{
    for (const auto& v : variants)
        WriteVariant(v);
}

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
