// Copyright (c) 2014, Pacific Biosciences of California, Inc.
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

#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include <string>
#include <algorithm>
#include <vector>

#include "FastaEntry.h"

namespace PBSeqAnalysis {
namespace PBChimera {

class FastaReader
{
public:  // structors
    // Default constructor
    FastaReader() = default;
    // Move constructor
    FastaReader(FastaReader&& src) = delete;
    // Copy constructor
    FastaReader(const FastaReader& src) = delete;
    // Move assignment constructor
    FastaReader& operator=(FastaReader&& rhs) = delete;
    // Copy assignment constructor
    FastaReader& operator=(const FastaReader& rhs) = delete;
    // Destructor
    ~FastaReader() = default;

public:  // non-modifying methods
    static std::shared_ptr<std::vector<FastaEntry>> ReadRecords(
            const char* filename)
    {
        auto output = std::shared_ptr<std::vector<FastaEntry>>(
                new std::vector<FastaEntry>);

        seqan::SeqFileIn seqFileIn;
        if (!seqan::open(seqFileIn, filename))
        {
            std::cout << "ERROR: Could not open the file.\n";
            return output;
        }

        seqan::StringSet<seqan::CharString> ids;
        seqan::StringSet<seqan::Dna5String> seqs;

        try
        {
            seqan::readRecords(ids, seqs, seqFileIn);
        }
        catch (seqan::Exception const & e)
        {
            std::cout << "ERROR: " << e.what() << std::endl;
            return output;
        }

        std::string id;
        seqan::Dna5String seq;
        uint32_t size;
        for (unsigned i = 0; i < seqan::length(ids); ++i)
        {
            id = std::string(seqan::toCString(ids[i]));
            seq = seqs[i];
            std::string numReadsStr = split(id, '_')[3];
            size = std::stoi(numReadsStr.substr(8));

            output->emplace_back(id, seq, seqan::length(seq), size);
        }

        std::sort(output->begin(), output->end(), FastaEntrySizeSort());

        return output;
    }

    static std::vector<std::string> split(const std::string &s, char delim) {
        std::vector<std::string> elems;
        std::stringstream ss(s);
        std::string item;
        while (std::getline(ss, item, delim)) {
            elems.push_back(item);
        }
        return elems;
    }
};

}  // namespace PBChimera
}  // namespace PBSeqAnalysis
