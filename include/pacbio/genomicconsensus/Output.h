// Copyright (c) 2017, Pacific Biosciences of California, Inc.
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

#pragma once

#include <iostream>
#include <string>
#include <vector>

#include <pbbam/FastaSequenceQuery.h>
#include <pbbam/QualityValues.h>
#include <pbcopper/utility/MoveAppend.h>

#include <pacbio/genomicconsensus/Settings.h>
#include <pacbio/genomicconsensus/WindowResult.h>

namespace PacBio {
namespace GenomicConsensus {

class Output
{
public:
    explicit Output(const Settings& settings);

    Output() = delete;
    Output(const Output&) = delete;
    Output(Output&&) = default;
    Output& operator=(const Output&) = delete;
    Output& operator=(Output&&) = default;
    ~Output() = default;

public:
    void AddResult(WindowResult result);

private:
    void MaybeFlushContig(const std::string& name);
    void RecordResult(WindowResult&& result);

private:
    Settings settings_;
    std::map<std::string, size_t> contigBasesRemaining_;
    std::map<std::string, std::vector<Consensus>> contigConsensi_;
    std::map<std::string, std::vector<Variant>> contigVariants_;
};

inline void PrintVcfLite(const Variant& v, std::ostream& out = std::cout)
{
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
            if (v.refSeq == v.readSeq1) {
                alt = v.readSeq2.get();
            } else if (v.refSeq == v.readSeq2.get()) {
                alt = v.readSeq1;
            } else {
                alt = v.readSeq1;
                alt.push_back(',');
                alt += v.readSeq2.get();
            }
        } else
            alt = v.readSeq1;
    }

    out << v.refName << '\t' << pos << "\t.\t" << ref << '\t' << alt << '\t'
        << v.confidence.get()  // force interpretation of 8-bit as num
        << "\tPASS\n";
}

inline Output::Output(const Settings& settings) : settings_{settings}
{
    // initialize maps
    PacBio::BAM::FastaSequenceQuery fasta{settings_.referenceFilename};
    for (const auto& seq : fasta) {
        contigBasesRemaining_[seq.Name()] = seq.Bases().size();
        contigConsensi_[seq.Name()] = std::vector<Consensus>{};
        contigVariants_[seq.Name()] = std::vector<Variant>{};
    }
}

inline void Output::AddResult(WindowResult result)
{
    const auto name = result.css.window.name;
    RecordResult(std::move(result));
    MaybeFlushContig(name);
}

inline void Output::MaybeFlushContig(const std::string& name)
{
    if (contigBasesRemaining_.at(name) == 0) {
        const auto& css = Consensus::Join(contigConsensi_.at(name));
        std::cout << "CSS:\n"
                  << css.sequence << '\n'
                  << "+\n"
                  << PacBio::BAM::QualityValues(css.confidence).Fastq() << '\n'
                  << "Variants:\n";
        for (const auto& v : contigVariants_.at(name))
            PrintVcfLite(v);
        std::cout << std::endl;

        contigConsensi_.at(name).clear();
        contigVariants_.at(name).clear();
    }
}

inline void Output::RecordResult(WindowResult&& result)
{
    const auto name = result.css.window.name;
    contigBasesRemaining_.at(name) -= result.css.window.Length();
    contigConsensi_.at(name).emplace_back(std::move(result.css));
    PacBio::Utility::MoveAppend(std::move(result.variants), contigVariants_.at(name));
}

}  // namespace GenomicConsensus
}  // namespace PacBio
