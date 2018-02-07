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

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <pacbio/genomicconsensus/NoCallStyle.h>
#include <pacbio/genomicconsensus/ReferenceWindow.h>

namespace PacBio {
namespace GenomicConsensus {

struct Consensus
{
    ReferenceWindow window;
    std::string sequence;
    std::vector<uint8_t> confidence;

    static Consensus NoCallConsensus(const NoCallStyle style, const ReferenceWindow& window,
                                     const std::string& refSeq);

    static Consensus Join(std::vector<Consensus> subconsensi);

    bool operator<(const Consensus& other) const { return window < other.window; }
};

inline bool AreContiguous(const std::vector<ReferenceWindow>& windows)
{
    //
    // Predicate that determines whether the reference/scaffold windows
    // are contiguous.
    //

    std::string lastName;
    size_t lastEnd = 0;
    for (const auto& win : windows) {
        if ((!lastName.empty() && win.name != lastName) ||
            (lastEnd != 0 && win.Start() != lastEnd)) {
            return false;
        }
        lastEnd = win.End();
        lastName = win.name;
    }
    return true;
}

inline Consensus Consensus::Join(std::vector<Consensus> subconsensi)
{
    //
    // [Consensus] -> Consensus
    //
    // String together all the consensus objects into a single consensus.
    // Will raise a ValueError if the reference windows are not
    // contiguous.
    //

    if (subconsensi.empty()) throw std::runtime_error("cannot join empty Consensus chunk list");
    std::sort(subconsensi.begin(), subconsensi.end());
    std::vector<ReferenceWindow> windows;
    for (const auto& c : subconsensi)
        windows.push_back(c.window);
    if (!AreContiguous(windows)) throw std::runtime_error("Consensus chunks must be contiguous");

    std::string joinedSeq;
    std::vector<uint8_t> joinedConfidence;
    for (const auto& c : subconsensi) {
        joinedSeq += c.sequence;
        joinedConfidence.insert(joinedConfidence.end(), c.confidence.begin(), c.confidence.end());
    }

    return Consensus{
        ReferenceWindow{subconsensi.front().window.name,
                        {subconsensi.front().window.Start(), subconsensi.back().window.End()}},
        std::move(joinedSeq), std::move(joinedConfidence)};
}

inline Consensus Consensus::NoCallConsensus(const NoCallStyle style, const ReferenceWindow& window,
                                            const std::string& refSeq)
{
    const auto length = refSeq.size();
    switch (style) {
        case (NoCallStyle::NO_CALL): {
            return Consensus{window, std::string(length, 'N'), std::vector<uint8_t>(length, 0)};
        }
        case (NoCallStyle::REFERENCE): {
            return Consensus{window, refSeq, std::vector<uint8_t>(length, 0)};
        }
        case (NoCallStyle::LOWERCASE_REFERENCE): {
            std::string seq;
            seq.reserve(refSeq.size());
            std::transform(refSeq.begin(), refSeq.end(), std::back_inserter(seq),
                           [](const char c) { return std::tolower(c); });
            return Consensus{window, std::move(seq), std::vector<uint8_t>(length, 0)};
        }
        default: {
            throw(std::string{"Unknown reference base call style!"});
        }
    }
}

}  // namespace GenomicConsensus
}  // namespace PacBio
