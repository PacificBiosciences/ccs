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

#include <pacbio/genomicconsensus/experimental/Output.h>

#include <iostream>

#include <pbcopper/logging/Logging.h>

#include <pacbio/genomicconsensus/experimental/Input.h>
#include <pacbio/genomicconsensus/experimental/Workflow.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

Output::Output(const Settings& settings) : settings_{settings}
{
    const Input input{settings_};
    const auto refWindows = input.ReferenceWindows(false);

    for (const auto& window : refWindows) {
        refWindows_[window.name] = window;
        processedBasesPerRef_[window.name] = 0;
        expectedBasesPerRef_[window.name] = window.Length();
        consensiPerRef_[window.name] = {};
        variantsPerRef_[window.name] = {};
    }

    // init writers from settings
    if (!settings_.fastaFilename.empty()) fasta_ = std::make_unique<FastaWriter>(settings_);
    if (!settings_.fastqFilename.empty()) fastq_ = std::make_unique<FastqWriter>(settings_);
    if (!settings_.gffFilename.empty()) gff_ = std::make_unique<GffWriter>(settings_, refWindows);
    if (!settings_.vcfFilename.empty()) vcf_ = std::make_unique<VcfWriter>(settings_, refWindows);
}

void Output::AddResult(WindowResult result)
{
    PBLOG_INFO << "Adding result for " << result.css.window;

    const auto window = result.css.window;
    consensiPerRef_[window.name].emplace_back(std::move(result.css));
    variantsPerRef_[window.name];  // append list
    processedBasesPerRef_[window.name] += window.Length();
    MaybeFlushContig(window.name);
}

void Output::MaybeFlushContig(const std::string& refName)
{
    const auto basesProcessed = processedBasesPerRef_[refName];
    const auto expectedBases = expectedBasesPerRef_[refName];
    if (basesProcessed == expectedBases) {
        // This contig is done, so we can dump to file & delete data
        if (gff_ || vcf_) {
            // sort & write variants
            auto& variants = variantsPerRef_[refName];
            std::sort(variants.begin(), variants.end());
            if (gff_) gff_->WriteVariants(variants);
            if (vcf_) vcf_->WriteVariants(variants);
        }

        // TODO: this doesn't actually free up memory
        variantsPerRef_[refName].clear();

        //
        // If the user asked to analyze a window or a set of
        // windows, we output a FAST[AQ] contig per analyzed
        // window.  Otherwise we output a fasta contig per
        // reference contig.
        //
        // We try to be intelligent about naming the output
        // contigs, to include window information where applicable.
        //
        for (const auto& window : Workflow::EnumerateWindows(refName, settings_)) {
            std::string recordName = window.name;
            if (window != refWindows_[window.name]) {
                recordName +=
                    '_' + std::to_string(window.Start()) + '_' + std::to_string(window.End());
            }

            std::string algoName;
            switch (settings_.mode) {
                case ConsensusMode::ARROW:
                    algoName = "arrow";
                    break;
                case ConsensusMode::PLURALITY:
                    algoName = "plurality";
                    break;
                case ConsensusMode::POA:
                    algoName = "poa";
                    break;
                default:
                    throw std::runtime_error("unknown consensus mode");
            }
            std::string cssName = recordName + '|' + algoName;

            std::vector<Consensus> consensiInThisWindow;
            for (const auto& consensus : consensiPerRef_[refName]) {
                if (Overlap(consensus.window, window)) consensiInThisWindow.push_back(consensus);
            }
            auto css = Consensus::Join(consensiInThisWindow);

            if (fasta_) fasta_->Write(cssName, css.sequence);
            if (fastq_) fastq_->Write(cssName, css.sequence, css.confidence);
        }

        // TODO: this doesn't actually free up memory
        consensiPerRef_[refName].clear();
    }
}

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
