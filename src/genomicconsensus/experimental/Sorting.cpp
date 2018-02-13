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

#include <pacbio/genomicconsensus/experimental/Sorting.h>

#include <algorithm>
#include <stdexcept>
#include <tuple>

#include <pacbio/genomicconsensus/experimental/ReferenceWindow.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

void SortReadsInWindow(std::vector<PacBio::BAM::BamRecord>* const reads,
                       const ReferenceWindow& window, const SortingStrategy strategy)
{
    using BamRecord = PacBio::BAM::BamRecord;

    auto lengthInWindow = [&window](const BamRecord& read) {
        const auto start = std::max(static_cast<size_t>(read.ReferenceStart()), window.Start());
        const auto end = std::min(static_cast<size_t>(read.ReferenceEnd()), window.End());
        return end - start;
    };

    switch (strategy) {
        case SortingStrategy::LONGEST_AND_STRAND_BALANCED: {
            // lexsort by read start, then end
            std::stable_sort(reads->begin(), reads->end(),
                             [](const BamRecord& lhs, const BamRecord& rhs) {
                                 const auto lhsStart = lhs.ReferenceStart();
                                 const auto rhsStart = rhs.ReferenceStart();
                                 const auto lhsEnd = lhs.ReferenceEnd();
                                 const auto rhsEnd = rhs.ReferenceEnd();
                                 return std::tie(lhsStart, lhsEnd) < std::tie(rhsStart, rhsEnd);
                             });

            // sort (descending) by reads' length in window
            std::stable_sort(reads->begin(), reads->end(),
                             [&window, lengthInWindow](const BamRecord& lhs, const BamRecord& rhs) {
                                 const auto lhsLength = lengthInWindow(lhs);
                                 const auto rhsLength = lengthInWindow(rhs);
                                 return lhsLength > rhsLength;
                             });
            break;
        }
        case SortingStrategy::LONGEST: {
            std::stable_sort(reads->begin(), reads->end(),
                             [&window, lengthInWindow](const PacBio::BAM::BamRecord& lhs,
                                                       const PacBio::BAM::BamRecord& rhs) {
                                 const auto lhsLength = lengthInWindow(lhs);
                                 const auto rhsLength = lengthInWindow(rhs);
                                 return rhsLength < lhsLength;
                             });
            break;
        }
        case SortingStrategy::SPANNING: {
            reads->erase(
                std::remove_if(reads->begin(), reads->end(),
                               [&window, lengthInWindow](const PacBio::BAM::BamRecord& record) {
                                   const auto length = lengthInWindow(record);
                                   return length != window.Length();
                               }),
                reads->end());
            break;
        }
        case SortingStrategy::FILE_ORDER: {
            // no sorting necessary
            break;
        }
        default:
            throw std::runtime_error("unexpected sorting strategy");
    }
}

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
