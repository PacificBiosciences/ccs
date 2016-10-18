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

#include <numeric>
#include <vector>

#include <pacbio/data/ArrayRead.h>
#include <pacbio/data/MSAColumn.h>

#include <pacbio/data/MSA.h>

namespace PacBio {
namespace Data {
MSA::MSA(const std::vector<Data::ArrayRead>& reads)
{
    // Determine left- and right-most positions
    for (const auto& r : reads) {
        beginPos = std::min(beginPos, r.ReferenceStart());
        endPos = std::max(endPos, r.ReferenceEnd());
    }
    // Fill counts
    FillCounts(reads);
}

void MSA::FillCounts(const std::vector<ArrayRead>& reads)
{
    // Prepare 2D array for counts
    counts.resize(endPos - beginPos);

    // Fill in counts
    for (const auto& r : reads) {
        int pos = r.ReferenceStart() - beginPos;
        assert(pos >= 0);
        std::string insertion;
        auto CheckInsertion = [&insertion, this, &pos]() {
            if (insertion.empty()) return;
            counts[pos].insertions[insertion]++;
            insertion = "";
        };
        for (const auto& b : r.Bases) {
            switch (b.Cigar) {
                case 'X':
                case '=':
                    CheckInsertion();
                    counts[pos][b.Nucleotide]++;
                    ++pos;
                    break;
                case 'D':
                    CheckInsertion();
                    counts[pos]['N']++;
                    ++pos;
                    break;
                case 'I':
                    insertion += b.Nucleotide;
                    break;
                case 'P':
                    CheckInsertion();
                    break;
                default:
                    throw std::runtime_error("Unexpected cigar " + std::to_string(b.Cigar));
            }
        }
    }
}
}  // namespace Data
}  // namespace PacBio
