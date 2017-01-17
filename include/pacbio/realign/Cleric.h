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

#pragma once

#include <exception>
#include <fstream>
#include <string>
#include <vector>

#include <pbbam/BamReader.h>
#include <pbbam/BamRecord.h>
#include <pbbam/BamWriter.h>
#include <pbbam/FastaReader.h>

namespace PacBio {
namespace Realign {

class Cleric
{
public:
    Cleric(const std::string& alignmentPath, const std::string& outputFile,
           const std::string& fromReference, const std::string& fromReferenceName,
           const std::string& toReference, const std::string& toReferenceName)
        : alignmentPath_(alignmentPath)
        , fromReferenceName_(fromReferenceName)
        , toReferenceName_(toReferenceName)
    {
        Align(fromReference, toReference, &fromReferenceSequence_, &toReferenceSequence_);
        Convert(outputFile);
    }

private:
    void Convert(const std::string& outputFile);
    void Align(const std::string& fromReference, const std::string& toReference,
               std::string* fromReferenceAligned, std::string* toReferenceAligned);

private:
    //clang-format off
    const BAM::CigarOperation newMatch_ =
        BAM::CigarOperation(BAM::CigarOperationType::SEQUENCE_MATCH, 1);
    const BAM::CigarOperation newDel_ = BAM::CigarOperation(BAM::CigarOperationType::DELETION, 1);
    const BAM::CigarOperation newIns_ = BAM::CigarOperation(BAM::CigarOperationType::INSERTION, 1);
    const BAM::CigarOperation newPad_ = BAM::CigarOperation(BAM::CigarOperationType::PADDING, 1);
    const BAM::CigarOperation newSoft_ = BAM::CigarOperation(BAM::CigarOperationType::SOFT_CLIP, 1);
    const BAM::CigarOperation newHard_ = BAM::CigarOperation(BAM::CigarOperationType::HARD_CLIP, 1);
    //clang-format on

private:
    const std::string alignmentPath_;
    std::string fromReferenceSequence_;
    std::string fromReferenceName_;
    std::string toReferenceSequence_;
    std::string toReferenceName_;

    std::string toReferenceGapless_;
    std::map<int, int> fasta_pos_to_sam_pos;

    std::string fromReferenceGapless_;
    std::map<int, int> sam_pos_to_fasta_pos;
};
}
}  // ::PacBio::Realign