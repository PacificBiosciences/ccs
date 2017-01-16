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

#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <pbcopper/cli/CLI.h>

#include <pbbam/BamReader.h>
#include <pbbam/FastaReader.h>

#include <pacbio/io/Utility.h>
#include <pacbio/realign/Cleric.h>
#include <pacbio/realign/ClericSettings.h>

namespace {
using namespace PacBio::Realign;

static void ParseInputFiles(const std::vector<std::string>& inputs, std::string* bamPath,
                            std::string* fromReference, std::string* fromReferenceName,
                            std::string* toReference, std::string* toReferenceName)
{
    using namespace PacBio::BAM;
    std::vector<std::string> fastaPaths;
    for (const auto& i : inputs) {
        try {
            BamReader reader(i);
            if (!bamPath->empty()) throw std::runtime_error("Only one BAM input is allowed!");
            *bamPath = i;
            if (reader.Header().Sequences().empty())
                throw std::runtime_error("Could not find reference sequence name");
            *fromReferenceName = reader.Header().Sequences().begin()->Name();
        } catch (...) {
            // If this is trigerred, the input file is not a BAM file.
            fastaPaths.push_back(i);
        }
    }

    for (const auto& fasta : fastaPaths) {
        FastaReader msaReader(fasta);

        FastaSequence f;
        while (msaReader.GetNext(f)) {
            if (f.Name() == *fromReferenceName) {
                if (fromReference->empty()) {
                    *fromReference = f.Bases();
                    std::transform(fromReference->begin(), fromReference->end(),
                                   fromReference->begin(), ::toupper);
                } else
                    throw std::runtime_error("Multiple original references provided!");
            } else if (toReference->empty()) {
                *toReference = f.Bases();
                std::transform(toReference->begin(), toReference->end(), toReference->begin(),
                               ::toupper);
                *toReferenceName = f.Name();
            } else {
                throw std::runtime_error("Multiple target references provided!");
            }
        }
    }
}

static int Runner(const PacBio::CLI::Results& options)
{
    // Check args size, as pbcopper does not enforce the correct number
    if (options.PositionalArguments().empty()) {
        std::cerr << "ERROR: Please provide BAM input, see --help" << std::endl;
        return EXIT_FAILURE;
    }
    if (options.PositionalArguments().size() >= 4) {
        std::cerr
            << "ERROR: Please provide _one_ BAM input and maximal _two_ FASTA files, see --help"
            << std::endl;
        return EXIT_FAILURE;
    }

    // Parse options
    ClericSettings settings(options);

    std::string bamPath;
    std::string fromReference;
    std::string fromReferenceName;
    std::string toReference;
    std::string toReferenceName;
    ParseInputFiles(settings.InputFiles, &bamPath, &fromReference, &fromReferenceName, &toReference,
                    &toReferenceName);

    std::string output;
    if (settings.OutputPrefix.empty())
        output = PacBio::IO::FilePrefix(bamPath) + "_cleric.bam";
    else
        output = settings.OutputPrefix + ".bam";

    Cleric cleric(bamPath, output, fromReference, fromReferenceName, toReference, toReferenceName);

    return EXIT_SUCCESS;
}
};

// Entry point
int main(int argc, char* argv[])
{
    return PacBio::CLI::Run(argc, argv, ClericSettings::CreateCLI(), &Runner);
}