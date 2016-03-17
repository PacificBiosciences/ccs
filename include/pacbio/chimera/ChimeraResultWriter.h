// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
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

// Author: Brett Bowman

#pragma once

#include <vector>

#include "ChimeraLabel.h"

namespace PacBio {
namespace Chimera {

///
/// \brief Write out to file a series ChimeraLabels in human-readable form
///
class ChimeraResultWriter
{
public:
    ChimeraResultWriter(const std::string &csvFile) : csvWriter_(csvFile)
    {
        ValidateCsvFile(csvFile);
        csvWriter_ << "SequenceId,IsChimera,ChimeraScore,"
                   << "ParentSequenceA,ParentSequenceB,CrossoverPosition" << std::endl;
    }

    void WriteResult(const ChimeraLabel &label) { csvWriter_ << label << std::endl; }
    void WriteResults(const std::vector<ChimeraLabel> &labels)
    {
        for (const auto &label : labels)
            WriteResult(label);
    }

private:
    static void ValidateCsvFile(const std::string &filename)
    {
        const std::string ext = filename.substr(filename.find_last_of(".") + 1);

        if (ext != "csv")
            throw std::invalid_argument("invalid sequence file type, only CSVs supported");
    }

private:
    std::ofstream csvWriter_;
};

}  // namespace Chimera
}  // namespace PacBio
