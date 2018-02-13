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

#include <pacbio/genomicconsensus/experimental/arrow/ArrowModel.h>

#include <pbbam/BamRecord.h>

#include <pacbio/genomicconsensus/experimental/arrow/Arrow.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

WindowResult ArrowModel::ConsensusAndVariantsFromWindow(
    const Input& /*input*/, const std::vector<PacBio::BAM::BamRecord>& reads,
    const ReferenceWindow& window, const std::string& refSeq, const Settings& settings) const
{
    Consensus css;
    boost::optional<std::vector<uint8_t>> effectiveSiteCoverage;
    if (settings.reportEffectiveCoverage) {
        std::vector<PacBio::BAM::BamRecord> readsUsed;
        css = Arrow::ConsensusForAlignments(window, refSeq, reads, settings, &readsUsed);
        effectiveSiteCoverage = Arrow::CoverageInWindow(window, readsUsed);
    } else
        css = Arrow::ConsensusForAlignments(window, refSeq, reads, settings, nullptr);

    const auto siteCoverage = Arrow::CoverageInWindow(window, reads);
    return Arrow::VariantsFromConsensus(window, refSeq, css, siteCoverage, effectiveSiteCoverage,
                                        settings);
}

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
