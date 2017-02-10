// Copyright (c) 2016, Pacific Biosciences of California, Inc.
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

// Inspired by work of David Seifert

#include <pacbio/align/SimdAlignment.h>

#include <pacbio/realign/Cleric.h>

namespace PacBio {
namespace Realign {
void Cleric::Align(const std::string& fromReference, const std::string& toReference,
                   std::string* fromReferenceAligned, std::string* toReferenceAligned)
{
    auto align = Align::SimdNeedleWunschAlignment(fromReference, toReference);

    *fromReferenceAligned = align.Target;
    *toReferenceAligned = align.Query;
}

void Cleric::Convert(const std::string& outputFile)
{
    using namespace PacBio::BAM;

    BamReader in(alignmentPath_);

    if (in.Header().Sequences().empty())
        throw std::runtime_error("Could not find reference sequence name");

    const auto localFromReferenceName = in.Header().Sequences().begin()->Name();
    if (localFromReferenceName != fromReferenceName_)
        throw std::runtime_error("Internal error. Reference name mismatches");

    const auto RemoveGaps = [](const std::string& input) {
        std::string seq = input;
        seq.erase(std::remove(seq.begin(), seq.end(), '-'), seq.end());
        return seq;
    };
    toReferenceGapless_ = RemoveGaps(toReferenceSequence_);
    fromReferenceGapless_ = RemoveGaps(fromReferenceSequence_);

    const auto GetGaplessMapping = [](const std::string& refBases, std::map<int, int>* map) {
        int pos = 0;
        for (size_t i = 0; i < refBases.size(); ++i) {
            if (refBases.at(i) != '-') {
                if (map->find(i) == map->cend()) {
                    map->emplace(pos, i);
                }
                ++pos;
            }
        }
    };

    GetGaplessMapping(fromReferenceSequence_, &sam_pos_to_fasta_pos);
    GetGaplessMapping(toReferenceSequence_, &fasta_pos_to_sam_pos);

    BamHeader h = in.Header().DeepCopy();
    h.ClearSequences();
    h.AddSequence(SequenceInfo(toReferenceName_, std::to_string(toReferenceGapless_.size())));

    BamWriter out(outputFile, h);
    BamRecord read;
    while (in.GetNext(read)) {
        std::string source_str = fromReferenceSequence_;
        std::string dest_str = toReferenceSequence_;

        // Expand RLE cigar to flat vector
        std::string expanded_cigar_ops;
        for (const auto& c : read.CigarData(false))
            for (size_t i = 0; i < c.Length(); ++i)
                expanded_cigar_ops += c.Char();
        expanded_cigar_ops += "YZ";

        CigarOperation old_cigar_state;  // UNKNOW_OP
        CigarOperation new_cigar_state;  // UNKNOW_OP

        bool found_start = false;
        int pos_in_read = 0;
        int pos_in_cigar = 0;
        int pos_in_source_ref = sam_pos_to_fasta_pos.at(read.ReferenceStart());

        Cigar new_cigar_tuple;

        int need_to_clip_left = 0;
        int need_to_clip_right = 0;

        int new_sam_start = 0;
        int pos_in_dest_ref = 0;

        while (pos_in_cigar < static_cast<int>(expanded_cigar_ops.size())) {
            char op = expanded_cigar_ops.at(pos_in_cigar);

            CigarOperation new_state;  // UNKNOWN_OP

            bool isFirstCigarAfterEnd = false;
            bool isSecondCigarAfterEnd = false;

            switch (op) {
                case 'M':
                case '=':
                case 'X':
                    if (!found_start) {
                        if (source_str.at(pos_in_source_ref) == '-') {
                            // Dest:   A---AAA
                            // Source: AAA-AAA
                            // Read:     A-AAA
                            //            ^

                            ++pos_in_source_ref;
                            continue;
                        }

                        // don't have a start POS yet
                        if (fasta_pos_to_sam_pos.find(pos_in_source_ref) !=
                            fasta_pos_to_sam_pos.cend()) {
                            new_sam_start = fasta_pos_to_sam_pos.at(pos_in_source_ref);
                            // Dest:   ---AAA
                            // Source: AAAAAA
                            // Read:      AAA
                            //            ^

                            new_state = newMatch_;
                            pos_in_dest_ref = pos_in_source_ref;
                            found_start = true;
                            ++pos_in_dest_ref;
                        } else {
                            // Dest:   ----AA
                            // Source: AAAAAA
                            // Read:      AAA
                            //            ^

                            // left Clip
                            new_state = newSoft_;
                        }

                        ++pos_in_cigar;
                        ++pos_in_read;
                        ++pos_in_source_ref;
                    } else {
                        if (source_str.at(pos_in_source_ref) == '-') {
                            if (dest_str.at(pos_in_dest_ref) == '-') {
                                // Dest:   AAA-AAA
                                // Source: AAA-AAA
                                // Read:   AAA-AAA
                                //            ^

                                ++pos_in_source_ref;
                                ++pos_in_dest_ref;
                                continue;
                            } else {
                                // Dest:   AAAAAAA
                                // Source: AAA-AAA
                                // Read:   AAA-AAA
                                //            ^

                                // Deletion
                                new_state = newDel_;

                                ++pos_in_source_ref;
                                ++pos_in_dest_ref;
                            }
                        } else {
                            if (dest_str.at(pos_in_dest_ref) == '-') {
                                // Dest:   AAA-AAA
                                // Source: AAAAAAA
                                // Read:   AAAAAAA
                                //            ^

                                // Insertion
                                new_state = newIns_;

                                ++pos_in_source_ref;
                                ++pos_in_dest_ref;
                                ++pos_in_cigar;
                                ++pos_in_read;
                            } else {
                                // Dest:   AAAAAAA
                                // Source: AAAAAAA
                                // Read:   AAAAAAA
                                //            ^

                                new_state = newMatch_;

                                ++pos_in_source_ref;
                                ++pos_in_dest_ref;
                                ++pos_in_cigar;
                                ++pos_in_read;
                            }
                        }
                    }
                    break;
                case 'I':
                    if (!found_start) {
                        if (source_str.at(pos_in_source_ref) == '-') {
                            // Dest:   A---AAA
                            // Source: AAA-AAA
                            // Read:     AAAAA
                            //            ^

                            ++pos_in_source_ref;
                            continue;
                        }

                        // Dest:   -- AAA
                        // Source: AA AAA
                        // Read:    AGAAA
                        //           ^

                        // left Clip
                        new_state = newSoft_;

                        ++pos_in_cigar;
                        ++pos_in_read;
                    } else {
                        if (source_str.at(pos_in_source_ref) == '-') {
                            if (dest_str.at(pos_in_dest_ref) == '-') {
                                // Dest:   AAA -AAA
                                // Source: AAA -AAA
                                // Read:   AAAA AAA
                                //            ^

                                ++pos_in_source_ref;
                                ++pos_in_dest_ref;
                                continue;
                            } else {
                                // Dest:   AAA AAAA
                                // Source: AAA -AAA
                                // Read:   AAAA AAA
                                //            ^

                                new_state = newMatch_;

                                ++pos_in_source_ref;
                                ++pos_in_dest_ref;
                                ++pos_in_cigar;
                                ++pos_in_read;
                            }
                        } else {
                            // Dest:   AAA -AAA
                            // Source: AAA AAAA
                            // Read:   AAAA AAA
                            //            ^
                            // OR
                            // Dest:   AAA AAA
                            // Source: AAA AAA
                            // Read:   AAAAAAA
                            //            ^

                            // Insertion
                            new_state = newIns_;

                            ++pos_in_cigar;
                            ++pos_in_read;
                        }
                    }
                    break;
                case 'N':
                case 'D':
                    if (!found_start) {
                        if (source_str.at(pos_in_source_ref) == '-') {
                            // Dest:   A---AAA
                            // Source: AAA-AAA
                            // Read:     A--AA
                            //            ^

                            ++pos_in_source_ref;
                            continue;
                        }

                        // Dest:   ---AAA
                        // Source: AAAAAA
                        // Read:    A-AAA
                        //           ^

                        ++pos_in_cigar;
                        ++pos_in_source_ref;
                        continue;
                    } else {
                        // have start POS
                        if (source_str.at(pos_in_source_ref) == '-') {
                            if (dest_str.at(pos_in_dest_ref) == '-') {
                                // Dest:   AAA-AAA
                                // Source: AAA-AAA
                                // Read:   AAA-AAA
                                //            ^

                                ++pos_in_source_ref;
                                ++pos_in_dest_ref;
                                continue;
                            } else {
                                // Dest:   AAAAAAA
                                // Source: AAA-AAA
                                // Read:   AAA--AA
                                //            ^

                                // Deletion
                                new_state = newDel_;

                                ++pos_in_source_ref;
                                ++pos_in_dest_ref;
                            }
                        } else {
                            if (dest_str.at(pos_in_dest_ref) == '-') {
                                // Dest:   AAA-AAA
                                // Source: AAAAAAA
                                // Read:   AAA-AAA
                                //            ^

                                // Padded deletion
                                ++pos_in_source_ref;
                                ++pos_in_dest_ref;
                                ++pos_in_cigar;

                                new_state = newPad_;
                            } else {
                                // Dest:   AAAAAAA
                                // Source: AAAAAAA
                                // Read:   AAA-AAA
                                //            ^

                                // Deletion
                                new_state = newDel_;

                                ++pos_in_source_ref;
                                ++pos_in_dest_ref;
                                ++pos_in_cigar;
                            }
                        }
                    }
                    break;
                case 'S':
                    new_state = newSoft_;

                    ++pos_in_cigar;
                    ++pos_in_read;
                    break;
                case 'H':
                    new_state = CigarOperation(CigarOperationType::HARD_CLIP, 1);

                    ++pos_in_cigar;
                    break;
                case 'P':
                    if (found_start) {
                        // Dest:   ---AAA
                        // Source: AAAAAA
                        // Read:    A-AAA
                        //           ^

                        ++pos_in_cigar;
                        ++pos_in_source_ref;
                        continue;

                    } else {
                        // have start POS
                        if (source_str.at(pos_in_source_ref) == '-') {
                            if (dest_str.at(pos_in_dest_ref) == '-') {
                                // Dest:   AAA-AAA
                                // Source: AAA-AAA
                                // Read:   AAA-AAA
                                //            ^

                                // Padded deletion
                                ++pos_in_cigar;

                                new_state = newPad_;
                            } else {
                                // Dest:   AAAAAAA
                                // Source: AAA-AAA
                                // Read:   AAA--AA
                                //            ^

                                // Deletion
                                new_state = newDel_;

                                ++pos_in_cigar;
                                ++pos_in_source_ref;
                                ++pos_in_dest_ref;
                            }
                        } else {
                            // Dest:   AAA--AAA
                            // Source: AAAAAAAA
                            // Read:   AAA-AAAA
                            //            ^
                            // OR
                            // Dest:   AAA AAAA
                            // Source: AAA AAAA
                            // Read:   AAA-AAAA
                            //            ^

                            // Padded deletion
                            ++pos_in_cigar;

                            new_state = newPad_;
                        }
                    }
                    break;
                case 'Y':
                    ++pos_in_cigar;
                    isFirstCigarAfterEnd = true;
                    break;
                case 'Z':
                    ++pos_in_cigar;
                    isSecondCigarAfterEnd = true;
                    break;
                default:
                    throw std::runtime_error("UNKNOWN CIGAR");
            }

            // If we reached Z, we have processed the CIGAR and can push the
            // lastest cigar operation.
            if (isSecondCigarAfterEnd) new_cigar_tuple.push_back(old_cigar_state);

            if (new_state.Type() != new_cigar_state.Type()) {
                // I ...... Y (end)
                if (new_state.Type() == CigarOperationType::UNKNOWN_OP && isFirstCigarAfterEnd &&
                    new_cigar_state.Type() == CigarOperationType::INSERTION) {
                    new_cigar_state.Type(CigarOperationType::SOFT_CLIP);
                }

                // have to rewrite CIGAR tuples if (a D and I operations are adjacen) {
                // D + I
                if (old_cigar_state.Type() == CigarOperationType::DELETION &&
                    new_cigar_state.Type() == CigarOperationType::INSERTION) {
                    const int num_del = old_cigar_state.Length();
                    const int num_insert = new_cigar_state.Length();
                    const int num_match = std::min(num_del, num_insert);

                    if (num_del == num_insert) {
                        // Dest:   GC AA-- TC      GC AA TC
                        // Read:   GC --AA TC  ->  GC AA TC
                        //            DDII            MM
                        old_cigar_state = CigarOperation();
                        new_cigar_state =
                            CigarOperation(CigarOperationType::SEQUENCE_MATCH, num_match);

                    } else if (num_del > num_insert) {
                        // Dest:   GC AAA-- TC      GC AAA TC
                        // Read:   GC ---AA TC  ->  GC -AA TC
                        //            DDDII            DMM
                        old_cigar_state =
                            CigarOperation(CigarOperationType::DELETION, num_del - num_match);
                        new_cigar_state =
                            CigarOperation(CigarOperationType::SEQUENCE_MATCH, num_match);

                    } else {
                        // Dest:   GC AA--- TC      GC AA- TC
                        // Read:   GC --AAA TC  ->  GC AAA TC
                        //            DDIII            MMI
                        old_cigar_state =
                            CigarOperation(CigarOperationType::SEQUENCE_MATCH, num_match);
                        new_cigar_state =
                            CigarOperation(CigarOperationType::INSERTION, num_insert - num_match);
                    }
                }

                // I + D
                if (old_cigar_state.Type() == CigarOperationType::INSERTION &&
                    new_cigar_state.Type() == CigarOperationType::DELETION) {
                    const int num_insert = old_cigar_state.Length();
                    const int num_del = new_cigar_state.Length();
                    const int num_match = std::min(num_del, num_insert);

                    if (num_del == num_insert) {
                        // Dest:   GC --AA TC  ->  GC AA TC
                        // Read:   GC AA-- TC      GC AA TC
                        //            IIDD            MM
                        old_cigar_state = CigarOperation();
                        new_cigar_state =
                            CigarOperation(CigarOperationType::SEQUENCE_MATCH, num_match);

                    } else if (num_del > num_insert) {
                        // Dest:   GC --AAA TC  ->  GC AAA TC
                        // Read:   GC AA--- TC      GC AA- TC
                        //            IIDDD            MMD
                        old_cigar_state =
                            CigarOperation(CigarOperationType::SEQUENCE_MATCH, num_match);
                        new_cigar_state =
                            CigarOperation(CigarOperationType::DELETION, num_del - num_match);

                    } else {
                        // Dest:   GC ---AA TC  ->  GC -AA TC
                        // Read:   GC AAA-- TC      GC AAA TC
                        //            IIIDD            IMM
                        old_cigar_state =
                            CigarOperation(CigarOperationType::INSERTION, num_insert - num_match);
                        new_cigar_state =
                            CigarOperation(CigarOperationType::SEQUENCE_MATCH, num_match);
                    }
                }

                if ((old_cigar_state.Type() != CigarOperationType::UNKNOWN_OP)) {
                    new_cigar_tuple.push_back(old_cigar_state);
                }
                // swap old and new state
                old_cigar_state = new_cigar_state;
                new_cigar_state = CigarOperation(new_state.Type(), 1);
            } else {
                new_cigar_state.Length(new_cigar_state.Length() + 1);
            }
        }

        //////////////////////////////////////
        // POST-PROCESSING                  //
        //////////////////////////////////////
        // check left flanking region + merge M-M pairs
        int i = 0;
        while (i < static_cast<int>(new_cigar_tuple.size()) - 1) {
            CigarOperation left_op = new_cigar_tuple.at(i);
            CigarOperation right_op = new_cigar_tuple.at(i + 1);

            // clang-format off
            // M + M:
            if (left_op.Type() == CigarOperationType::SEQUENCE_MATCH && right_op.Type() == CigarOperationType::SEQUENCE_MATCH) {
                new_cigar_tuple[i] =
                    CigarOperation(CigarOperationType::SEQUENCE_MATCH, left_op.Length() + right_op.Length());
                new_cigar_tuple.erase(new_cigar_tuple.begin() + i + 1);
            }
            // S + I:
            else if (left_op.Type() == CigarOperationType::SOFT_CLIP && right_op.Type() == CigarOperationType::INSERTION) {
                new_cigar_tuple[i] = CigarOperation(
                    CigarOperationType::SOFT_CLIP, left_op.Length() + right_op.Length());
                new_cigar_tuple.erase(new_cigar_tuple.begin() + i + 1);
            }
            // S + D:
            else if (left_op.Type() == CigarOperationType::SOFT_CLIP && right_op.Type() == CigarOperationType::DELETION) {
                new_cigar_tuple[i] = CigarOperation(CigarOperationType::SOFT_CLIP, left_op.Length());
                new_cigar_tuple.erase(new_cigar_tuple.begin() + i + 1);
            }
            // S + P:
            else if (left_op.Type() == CigarOperationType::SOFT_CLIP && right_op.Type() == CigarOperationType::PADDING) {
                new_cigar_tuple[i] = CigarOperation(CigarOperationType::SOFT_CLIP, left_op.Length());
                new_cigar_tuple.erase(new_cigar_tuple.begin() + i + 1);
            }
            // H + I:
            else if (left_op.Type() == CigarOperationType::HARD_CLIP && right_op.Type() == CigarOperationType::INSERTION) {
                new_cigar_tuple[i + 1] = CigarOperation(CigarOperationType::SOFT_CLIP, right_op.Length());
                ++i;
            }
            // H + D:
            else if (left_op.Type() == CigarOperationType::HARD_CLIP && right_op.Type() == CigarOperationType::DELETION) {
                new_cigar_tuple[i] = CigarOperation(CigarOperationType::HARD_CLIP, left_op.Length());
                new_cigar_tuple.erase(new_cigar_tuple.begin() + i + 1);
            }
            // H + P:
            else if (left_op.Type() == CigarOperationType::HARD_CLIP && right_op.Type() == CigarOperationType::PADDING) {
                new_cigar_tuple[i] = CigarOperation(CigarOperationType::HARD_CLIP, left_op.Length());
                new_cigar_tuple.erase(new_cigar_tuple.begin() + i + 1);
                // H + S:
                // } else if (left_op.Type() == CigarOperationType::HARD_CLIP && right_op.Type() == CigarOperationType::SOFT_CLIP) {
                //    ++i;
            } else {
                ++i;
            }
            // clang-format on
        }

        // check right flanking region
        i = new_cigar_tuple.size() - 2;
        // cant_stop = True
        while (i >= 0) {
            // cant_stop = False

            CigarOperation left_op = new_cigar_tuple.at(i);
            CigarOperation right_op = new_cigar_tuple.at(i + 1);

            if (left_op.Type() == CigarOperationType::SEQUENCE_MATCH) {
                // reached a match state, hence everything
                // before will be compliant
                break;
            }
            // I + S:
            if (left_op.Type() == CigarOperationType::INSERTION &&
                right_op.Type() == CigarOperationType::SOFT_CLIP) {
                // cant_stop = True
                new_cigar_tuple[i] = CigarOperation(CigarOperationType::SOFT_CLIP,
                                                    left_op.Length() + right_op.Length());
                new_cigar_tuple.erase(new_cigar_tuple.begin() + i + 1);
            }
            // D + S:
            else if (left_op.Type() == CigarOperationType::DELETION &&
                     right_op.Type() == CigarOperationType::SOFT_CLIP) {
                // cant_stop = True
                new_cigar_tuple[i] =
                    CigarOperation(CigarOperationType::SOFT_CLIP, right_op.Length());
                new_cigar_tuple.erase(new_cigar_tuple.begin() + i + 1);
            }
            // P + S:
            else if (left_op.Type() == CigarOperationType::PADDING &&
                     right_op.Type() == CigarOperationType::SOFT_CLIP) {
                // cant_stop = True
                new_cigar_tuple[i] =
                    CigarOperation(CigarOperationType::SOFT_CLIP, right_op.Length());
                new_cigar_tuple.erase(new_cigar_tuple.begin() + i + 1);
            }
            // I + H:
            else if (left_op.Type() == CigarOperationType::INSERTION &&
                     right_op.Type() == CigarOperationType::HARD_CLIP) {
                // cant_stop = True
                new_cigar_tuple[i] =
                    CigarOperation(CigarOperationType::SOFT_CLIP, left_op.Length());
            }
            // D + H:
            else if (left_op.Type() == CigarOperationType::DELETION &&
                     right_op.Type() == CigarOperationType::HARD_CLIP) {
                // cant_stop = True
                new_cigar_tuple[i] =
                    CigarOperation(CigarOperationType::HARD_CLIP, right_op.Length());
                new_cigar_tuple.erase(new_cigar_tuple.begin() + i + 1);
            }
            // P + H:
            else if (left_op.Type() == CigarOperationType::PADDING &&
                     right_op.Type() == CigarOperationType::HARD_CLIP) {
                // cant_stop = True
                new_cigar_tuple[i] =
                    CigarOperation(CigarOperationType::HARD_CLIP, right_op.Length());
                new_cigar_tuple.erase(new_cigar_tuple.begin() + i + 1);
            }
            // S + H:
            // } else if (left_op.Type() == CigarOperationType::SOFT_CLIP && right_op.Type() == CigarOperationType::HARD_CLIP) {
            //    //cant_stop = True
            //    pass
            --i;
        }

        std::string new_seq = read.Sequence(Orientation::GENOMIC);

        // calculate edit distance (and possibly replace match states)
        pos_in_read = 0;
        pos_in_dest_ref = new_sam_start;
        int new_edit_distance = 0;
        Cigar replace_cigar_tuple;

        const auto match_state_det = [](char read_base, char genome_base) {
            if (read_base == genome_base)
                return CigarOperationType::SEQUENCE_MATCH;
            else
                return CigarOperationType::SEQUENCE_MISMATCH;
        };

        int tlen = 0;
        for (const auto& op : new_cigar_tuple) {
            const CigarOperationType cigar_op = op.Type();
            const int cigar_op_count = op.Length();

            if (cigar_op == CigarOperationType::SEQUENCE_MATCH) {
                auto old_state = match_state_det(new_seq.at(pos_in_read),
                                                 toReferenceGapless_.at(pos_in_dest_ref));
                int count = 1;
                for (int i = 1; i < cigar_op_count; ++i) {
                    const auto next_state = match_state_det(
                        new_seq.at(pos_in_read + i), toReferenceGapless_.at(pos_in_dest_ref + i));
                    if (old_state != next_state) {
                        if (old_state == CigarOperationType::SEQUENCE_MISMATCH)
                            new_edit_distance += count;
                        replace_cigar_tuple.emplace_back(old_state, count);
                        old_state = next_state;
                        count = 1;
                    } else {
                        ++count;
                    }
                }

                if (old_state == CigarOperationType::SEQUENCE_MISMATCH) new_edit_distance += count;
                replace_cigar_tuple.emplace_back(old_state, count);

                tlen += cigar_op_count;
                pos_in_read += cigar_op_count;
                pos_in_dest_ref += cigar_op_count;
            } else if (cigar_op == CigarOperationType::INSERTION) {
                new_edit_distance += cigar_op_count;
                replace_cigar_tuple.emplace_back(cigar_op, cigar_op_count);
                pos_in_read += cigar_op_count;
            } else if (cigar_op == CigarOperationType::DELETION) {
                new_edit_distance += cigar_op_count;
                replace_cigar_tuple.emplace_back(cigar_op, cigar_op_count);
                tlen += cigar_op_count;
                pos_in_dest_ref += cigar_op_count;
            } else if (cigar_op == CigarOperationType::SOFT_CLIP) {
                replace_cigar_tuple.emplace_back(cigar_op, cigar_op_count);
                pos_in_read += cigar_op_count;
            } else if (cigar_op == CigarOperationType::HARD_CLIP ||
                       cigar_op == CigarOperationType::PADDING) {
                replace_cigar_tuple.emplace_back(cigar_op, cigar_op_count);
            } else {
                throw std::runtime_error("STATE should not occur " +
                                         std::to_string(static_cast<int>(cigar_op)));
            }
        }
        read.Impl().CigarData(replace_cigar_tuple);
        read.Impl().Position(new_sam_start);
        out.Write(read);
    }
}
}
}  // ::PacBio::Realign