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
    BeginEnd(reads);
    // Fill counts
    FillCounts(reads);
}
MSA::MSA(const std::vector<Data::ArrayRead>& reads, const MSA& prior)
{
    BeginEnd(reads);
    // Fill counts
    FillCounts(reads, prior);
}

void MSA::BeginEnd(const std::vector<Data::ArrayRead>& reads)
{
    // Determine left- and right-most positions
    for (const auto& r : reads) {
        beginPos = std::min(beginPos, r.ReferenceStart());
        endPos = std::max(endPos, r.ReferenceEnd());
    }
}

void MSA::FillCounts(const std::vector<ArrayRead>& reads)
{
    // Prepare 2D array for counts
    counts.resize(endPos - beginPos);
    int pos = beginPos;
    for (auto& c : counts)
        c.refPos = ++pos;  // output is 1-based, input is 0-based

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

void MSA::FillCounts(const std::vector<ArrayRead>& reads, const MSA& prior)
{
    // Prepare 2D array for counts
    counts.resize(endPos - beginPos);
    {
        int pos = beginPos;
        for (auto& c : counts)
            c.refPos = pos++;
    }

    struct InDel
    {
        InDel(const MSAColumn& c)
            : refPos(c.refPos), deletion(c.mask.at(4)), insertions(c.SignificantInsertions())
        {
        }

        bool Hit() { return deletion || !insertions.empty(); }
        int refPos = -1;
        bool deletion = false;
        std::vector<std::string> insertions;
    };
    std::vector<InDel> indels;
    for (auto& column : prior)
        indels.emplace_back(column);

    for (const auto& id : indels)
        if (id.deletion) std::cerr << id.refPos << " " << id.deletion << std::endl;

    std::map<int, int> offsets;
    std::map<int, int> delMap;
    // Fill in counts
    for (const auto& r : reads) {
        int pos = r.ReferenceStart() - beginPos;
        assert(pos >= 0);

        int indelOffset = 0;
        std::string insertion;
        int deletion = 0;

        auto CheckInsertion = [&insertion, &indels, &indelOffset, this, &pos]() {
            if (insertion.empty()) return;
            if (pos < static_cast<int>(indels.size())) {
                const auto& insertions = indels.at(pos).insertions;
                if (std::find(insertions.cbegin(), insertions.cend(), insertion) !=
                    insertions.cend()) {
                    if (insertion.size() % 3 != 0) indelOffset += insertion.size();
                    std::cerr << "Found insertion " << insertion << " at positition " << pos
                              << std::endl;
                }
            } else
                std::cerr << pos << " is outside of " << indels.size() << " for insertion "
                          << insertion << std::endl;
            insertion = "";
        };

        auto CheckDeletion = [&deletion, &indelOffset]() {
            if (deletion != 0) {
                if (deletion % 3 != 0) indelOffset -= deletion;
                deletion = 0;
            }
        };

        for (size_t i = 0; i < r.Bases.size(); ++i) {
            const auto& b = r.Bases.at(i);
            bool hasFullCodonFollowing = (i + 2) < r.Bases.size();
            switch (b.Cigar) {
                case 'X':
                case '=':
                    // Check for HP deletion
                    if ((deletion != 0) && (i + 1) < r.Bases.size() &&
                        b.Nucleotide == r.Bases.at(i + 1).Nucleotide) {
                        deletion = 0;
                    } else {  // deletion is not in an HP
                        CheckDeletion();
                    }
                    if (!insertion.empty()) {
                        bool hp = true;
                        for (size_t j = 0; j < insertion.size(); ++j) {
                            if (insertion[j] != b.Nucleotide) {
                                hp = false;
                                break;
                            }
                        }
                        if (!hp) CheckInsertion();
                    }
                    ++pos;
                    break;
                case 'D':
                    if (indels.at(pos).deletion) {
                        delMap[beginPos + pos + 1]++;
                        ++deletion;
                    }
                    CheckInsertion();
                    ++pos;
                    break;
                case 'I':
                    insertion += b.Nucleotide;
                    CheckDeletion();
                    break;
                case 'P':
                    CheckDeletion();
                    CheckInsertion();
                    break;
                default:
                    throw std::runtime_error("Unexpected cigar " + std::to_string(b.Cigar));
            }
        }
        offsets[indelOffset]++;

        // int pos = r.ReferenceStart() - beginPos;
        // assert(pos >= 0);
        // std::string insertion;
        // auto CheckInsertion = [&insertion, this, &pos]() {
        //     if (insertion.empty()) return;
        //     counts[pos].insertions[insertion]++;
        //     insertion = "";
        // };
        // for (const auto& b : r.Bases) {
        //     switch (b.Cigar) {
        //         case 'X':
        //         case '=':
        //             CheckInsertion();
        //             counts[pos][b.Nucleotide]++;
        //             ++pos;
        //             break;
        //         case 'D':
        //             CheckInsertion();
        //             counts[pos]['N']++;
        //             ++pos;
        //             break;
        //         case 'I':
        //             insertion += b.Nucleotide;
        //             break;
        //         case 'P':
        //             CheckInsertion();
        //             break;
        //         default:
        //             throw std::runtime_error("Unexpected cigar " + std::to_string(b.Cigar));
        //     }
        // }
    }
    std::cerr << "del" << std::endl;
    for (const auto& pos_count : delMap)
        std::cerr << pos_count.first << " - " << pos_count.second << std::endl;
    std::cerr << "offsets" << std::endl;
    for (const auto& offset_count : offsets)
        std::cerr << offset_count.first << " - " << offset_count.second << std::endl;
}
}  // namespace Data
}  // namespace PacBio
