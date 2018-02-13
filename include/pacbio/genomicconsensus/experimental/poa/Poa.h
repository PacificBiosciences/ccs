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

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include <pacbio/align/AffineAlignment.h>
#include <pacbio/align/PairwiseAlignment.h>
#include <pacbio/data/Interval.h>
#include <pacbio/genomicconsensus/experimental/Consensus.h>
#include <pacbio/genomicconsensus/experimental/Input.h>
#include <pacbio/genomicconsensus/experimental/Settings.h>
#include <pacbio/genomicconsensus/experimental/Utils.h>
#include <pacbio/genomicconsensus/experimental/WindowResult.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

struct Poa
{
    struct ConfidenceAndVariantResult
    {
        std::vector<uint8_t> confidence;
        std::vector<Variant> variants;
    };

    // refWindow, refSequence, poaCss, poaConfig.aligner
    static ConfidenceAndVariantResult ConfidenceAndVariants(const ReferenceWindow& window,
                                                            const std::string& refSeq,
                                                            const std::string& poaCss,
                                                            const Settings& settings)
    {
        //
        // Compute the confidence for each position, and compare
        // the consensus and reference in this window, returning a list of variants
        //

        using PacBio::Align::AlignAffineIupac;
        using PairwiseAlignment = PacBio::Align::PairwiseAlignment;

        std::unique_ptr<PairwiseAlignment> ga = nullptr;
        // TODO: switch on CLI-requested aligner, for now just go with affine (default)
        ga.reset(AlignAffineIupac(refSeq, poaCss));
        //        if (settings.aligner == "affine")
        //            ga.reset(AlignAffineIupac(refSeq, poaCss));
        //        else
        //            ga.reset(Align(refSeq, poaCss));
        assert(ga);

        auto confidence = std::vector<uint8_t>(poaCss.size(), 20);
        auto variants =
            VariantsFromAlignment(ga, window, confidence, std::vector<uint8_t>(), boost::none);
        return ConfidenceAndVariantResult{std::move(confidence), std::move(variants)};
    }

    static WindowResult ConsensusAndVariantsForAlignments(
        const ReferenceWindow& window, const std::string& refSeq,
        const std::vector<PacBio::BAM::BamRecord>& reads, const Settings& settings)
    {
        //
        // Call consensus on this interval---without subdividing the interval
        // further.
        //
        // Testable!
        //
        // Clipping has already been done!
        //

        const auto refStart = window.Start();
        const auto refEnd = window.End();

        // Compute the POA consensus, which is our initial guess, and
        // should typically be > 99.5% accurate
        auto fwdSequences = FilteredForwardSequences(reads, window);
        std::string poaCss;
        try {
            const auto p = MakePoaConsensus(std::move(fwdSequences), settings);
            poaCss = p->Sequence;
        } catch (std::exception&) {
            //TODO: log
            auto css =
                Consensus::NoCallConsensus(NoCallStyle::LOWERCASE_REFERENCE,  // double-check this
                                           window, refSeq);
            return WindowResult{std::move(css), std::vector<Variant>()};
        }

        auto confAndVars = ConfidenceAndVariants(window, refSeq, poaCss, settings);
        return WindowResult{Consensus{window, std::move(poaCss), std::move(confAndVars.confidence)},
                            std::move(confAndVars.variants)};
    }

    static WindowResult ConsensusAndVariantsForWindow(const Input& input,
                                                      const ReferenceWindow& window,
                                                      const std::string refSeq,
                                                      const Settings& settings)
    {
        //
        // High-level routine for calling the consensus for a
        // window of the genome given an alignment.
        //
        // Identifies the coverage contours of the window in order to
        // identify subintervals where a good consensus can be called.
        // Creates the desired "no evidence consensus" where there is
        // inadequate coverage.
        //

        using Interval = PacBio::Data::Interval;

        // TODO: fancy chunking

        std::vector<Consensus> subconsensi;
        std::vector<Variant> variants;

        std::vector<Interval> intervals{window.interval};
        for (const auto& interval : intervals) {

            const auto intervalRefSeq = refSeq.substr(interval.Left(), interval.Length());
            const ReferenceWindow subWindow{window.name, interval};

            auto reads = input.ReadsInWindow(subWindow);
            ClipReadsToWindow(&reads, subWindow);
            FilterAlignments(&reads, settings);

            // if enough coverage
            const size_t numSpanning = std::count_if(
                reads.begin(), reads.end(), [&interval](const PacBio::BAM::BamRecord& read) {
                    const auto readStart = static_cast<size_t>(read.ReferenceStart());
                    const auto readEnd = static_cast<size_t>(read.ReferenceEnd());
                    return readStart <= interval.Left() && interval.Right() <= readEnd;
                });

            Consensus css;
            if (numSpanning >= settings.minPoaCoverage) {

                auto tempWindowResult =
                    ConsensusAndVariantsForAlignments(subWindow, intervalRefSeq, reads, settings);

                // store window css
                css = std::move(tempWindowResult.css);

                // store window variants
                auto filteredVariants = FilterVariants(tempWindowResult.variants, settings);
                if (settings.annotateGFF) AnnotateVariants(&filteredVariants, reads);
                variants.insert(variants.end(), filteredVariants.begin(), filteredVariants.end());

                // maybe dump evidence

            } else  // not enough coverage
            {
                css = Consensus::NoCallConsensus(NoCallStyle::LOWERCASE_REFERENCE, subWindow,
                                                 intervalRefSeq);
            }

            // add result to our final consensi
            subconsensi.push_back(css);
        }

        return WindowResult{Consensus::Join(std::move(subconsensi)), std::move(variants)};
    }
};

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
