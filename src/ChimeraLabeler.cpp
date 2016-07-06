// Copyright (c) 2015, Pacific Biosciences of California, Inc.
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

#include <seqan/seq_io.h>

#include "pacbio/chimera/ChimeraLabeler.h"
#include "pacbio/chimera/ChimeraResultWriter.h"

using namespace PacBio::Chimera;

int main(int argc, char const **argv)
{
    using namespace seqan;

    if (argc != 2) return 1;  // Invalid number of arguments.

    // Get the input file
    std::string inputFile(argv[1]);
    SeqFileIn inputHandle(inputFile.c_str());

    // Parse the records
    StringSet<CharString> ids;
    StringSet<Dna5String> seqs;
    readRecords(ids, seqs, inputHandle);

    // Declare the vectors we'll use to actually perform the Chimera-labeling
    std::vector<std::string> idList;
    std::vector<Dna5String> seqList;

    for (size_t i = 0; i < length(ids); ++i) {
        idList.push_back(toCString(static_cast<CharString>(ids[i])));
        seqList.push_back(static_cast<Dna5String>(seqs[i]));
    }

    // Label the Records
    ChimeraLabeler chimeraLabeler(1.0f, 100, true);
    auto labels = chimeraLabeler.LabelChimeras(idList, seqList);

    // Display the results
    ChimeraResultWriter csvWriter("temp.csv");
    csvWriter.WriteResults(labels);

    return 0;
}
