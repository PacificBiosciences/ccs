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

// Author: David Seifert, Brett Bowman

#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <pbbam/../../src/SequenceUtils.h>
#include <pbbam/BamReader.h>
#include <pbbam/BamRecord.h>
#include <pbbam/BamRecordBuilder.h>
#include <pbbam/BamWriter.h>
#include <pbbam/Cigar.h>
#include <pbbam/FastaReader.h>
#include <pbbam/MD5.h>
#include <pbbam/ReadGroupInfo.h>
#include <pbbam/TagCollection.h>

#include <pacbio/Version.h>
#include <pacbio/consensus/ModelConfig.h>
#include "../ModelFactory.h"
#include "../Simulator.h"

using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace PacBio::Consensus;
using namespace PacBio::Data;

// these strings are part of the BAM header, they CANNOT contain newlines
const static std::string DESCRIPTION("Simulate genomic (sub)reads from an aligned PacBio BAM.");
const static std::string APPNAME("genomic_sim");

namespace {  // anonymous
inline PacBio::BAM::Cigar ConvertStatePathToCigar(const std::vector<MoveType>& statePath,
                                                  const std::string& tpl, const std::string& read)
{
    using namespace PacBio::BAM;
    Cigar result;

    bool foundStart = false;
    int32_t posOnTpl = 0, posOnRead = 0, editDistance = 0;

    CigarOperationType newType{CigarOperationType::UNKNOWN_OP};
    CigarOperationType curType{CigarOperationType::UNKNOWN_OP};
    uint32_t curLength = 0;

    for (const auto& state : statePath) {
        switch (state) {
            case MoveType::MATCH:
                if (!foundStart) {
                    foundStart = true;
                }

                if (tpl[posOnTpl] == read[posOnRead]) {
                    // match
                    newType = CigarOperationType::SEQUENCE_MATCH;
                } else {
                    // mismatch
                    newType = CigarOperationType::SEQUENCE_MISMATCH;
                    ++editDistance;
                }

                ++posOnTpl;
                ++posOnRead;
                break;

            case MoveType::BRANCH:
            case MoveType::STICK:
                if (foundStart) {
                    // encountered at least one match before
                    newType = CigarOperationType::INSERTION;
                    ++editDistance;
                } else {
                    // still in a clip region
                    newType = CigarOperationType::SOFT_CLIP;
                }
                ++posOnRead;
                break;

            case MoveType::DELETION:
                ++posOnTpl;
                if (foundStart) {
                    // encountered one match before
                    newType = CigarOperationType::DELETION;
                    ++editDistance;
                } else {
                    continue;
                }
                break;
        }

        if (curType != newType) {
            if (curType != CigarOperationType::UNKNOWN_OP) {
                result.emplace_back(curType, curLength);
            }
            curType = newType;
            curLength = 0;
        }
        ++curLength;
    }

    result.emplace_back(curType, curLength);

    return result;
}
}  // anonymous namespace

static BamHeader PrepareHeader(const std::string& cmdLine,
                               const std::vector<ReadGroupInfo>& inputReadgroups)
{
    ProgramInfo program{APPNAME + "-" + PacBio::UnanimityVersion()};
    program.Name(APPNAME)
        .CommandLine(APPNAME + " " + cmdLine)
        .Description(DESCRIPTION)
        .Version(PacBio::UnanimityVersion());

    std::random_device rd;
    std::default_random_engine rng{rd()};

    if (inputReadgroups.size() > 1) {
        std::cerr << APPNAME << " can currently only handle one readgroup per input bam\n";
        exit(EXIT_FAILURE);
    }

    std::ostringstream movieName;
    movieName << 'm' << std::uniform_int_distribution<uint32_t>{10000, 99999}(rng) << '_'
              << std::uniform_int_distribution<uint32_t>{100000, 999999}(rng) << '_'
              << std::uniform_int_distribution<uint32_t>{100000, 999999}(rng);

    // reuse input readgroups for platform information
    ReadGroupInfo newRg{inputReadgroups.front()};
    newRg.MovieName(movieName.str())
        .ReadType("SUBREAD")
        .PlatformModel(PlatformModelType::SEQUEL)
        .IpdCodec(FrameCodec::V1)
        .PulseWidthCodec(FrameCodec::V1)
        .Id(movieName.str(), "SUBREAD");

    BamHeader header;
    header.PacBioBamVersion("3.0.1")
        .SortOrder("unknown")
        .Version("1.5")
        .AddProgram(program)
        .ReadGroups(std::vector<ReadGroupInfo>{newRg});

    return header;
}

static int SimulateGenomicReads(const std::string& referenceFilename,
                                const std::string& inputFilename, const std::string& outputFilename,
                                unsigned int seed = 42)
{
    std::vector<FastaSequence> references = FastaReader::ReadAll(referenceFilename);

    BamReader reader{inputFilename};
    BamRecord inputRecord;

    BamHeader newHeader{PrepareHeader("genomic_sim", reader.Header().ReadGroups())};
    assert(newHeader.ReadGroups().size() == 1);
    const std::string newRg{newHeader.ReadGroups().front().Id()};

    // Create RC'd refs and add our references to the new header
    std::vector<FastaSequence> rcReferences;
    std::vector<size_t> refSizes;
    for (const auto& ref : references) {
        refSizes.push_back(ref.Bases().size());
        rcReferences.push_back(FastaSequence{ref.Name(), ReverseComplemented(ref.Bases())});
        SequenceInfo si{ref.Name(), std::to_string(ref.Bases().size())};
        si.Checksum(MD5Hash(ref.Bases()));
        newHeader.AddSequence(si);
    }

    std::vector<BamRecord> newRecords;
    std::default_random_engine rng{seed};

    int32_t zmw = 0;
    while (reader.GetNext(inputRecord)) {
        ++zmw;

        // 1. Check that our reference is still valid
        assert(static_cast<size_t>(inputRecord.ReferenceId()) < references.size());
        assert(static_cast<size_t>(inputRecord.ReferenceEnd()) <=
               references[inputRecord.ReferenceId()].Bases().size());

        // 2. extract the region to be simulated
        const int32_t refSpan = inputRecord.ReferenceEnd() - inputRecord.ReferenceStart();
        const bool isRevStrand = inputRecord.AlignedStrand() == Strand::REVERSE;
        std::string referenceSeq;
        if (isRevStrand) {
            const size_t refStart =
                refSizes[inputRecord.ReferenceId()] - inputRecord.ReferenceEnd();
            referenceSeq =
                rcReferences[inputRecord.ReferenceId()].Bases().substr(refStart, refSpan);
        } else {
            referenceSeq = references[inputRecord.ReferenceId()].Bases().substr(
                inputRecord.ReferenceStart(), refSpan);
        }

        // 3. simulate the new read
        std::unique_ptr<ModelConfig> currentModel = ModelFactory::Create(
            inputRecord.ReadGroup().SequencingChemistry(), inputRecord.SignalToNoise());
        std::pair<Read, std::vector<MoveType>> rawRead =
            currentModel->SimulateRead(&rng, referenceSeq, "");

        // 4. prepare new subread
        BamRecord newRecord{newHeader};

        // 5. orient the sequence, cigar, ipd and pulse-width data
        std::string newSeq = rawRead.first.Seq;
        Cigar newCigar{ConvertStatePathToCigar(rawRead.second, referenceSeq, newSeq)};
        std::vector<uint8_t> ipd(rawRead.first.IPD);
        std::vector<uint8_t> pw(rawRead.first.PulseWidth);
        if (isRevStrand) {
            ReverseComplement(newSeq);
            std::reverse(newCigar.begin(), newCigar.end());
            std::reverse(ipd.begin(), ipd.end());
            std::reverse(pw.begin(), pw.end());
        }

        // 5. fill out the read
        newRecord.ReadGroup(newRg)
            .IPD(Frames::Decode(ipd), FrameEncodingType::LOSSY)
            .NumPasses(1)
            .PulseWidth(Frames::Decode(pw), FrameEncodingType::LOSSY)
            .QueryStart(0)
            .QueryEnd(newSeq.length())
            .ReadAccuracy(0.8)
            .SignalToNoise(rawRead.first.SignalToNoise)
            .HoleNumber(zmw)
            .UpdateName();

        newRecord.Impl()
            .CigarData(newCigar)
            .Bin(0)
            .InsertSize(0)
            .MapQuality(254)
            .MatePosition(-1)
            .MateReferenceId(-1)
            .Position(inputRecord.ReferenceStart())
            .ReferenceId(inputRecord.ReferenceId())
            .SetMapped(true)
            .SetReverseStrand(isRevStrand)
            .SetSequenceAndQualities(newSeq);

        //  6. append the original read-name for record-keeping
        newRecord.Impl().AddTag("fn", inputRecord.FullName());

        newRecords.push_back(newRecord);
    }

    // write bam file
    BamWriter newWriter{outputFilename, newHeader, BamWriter::BestCompression};
    for (const auto& i : newRecords) {
        newWriter.Write(i);
    }

    return EXIT_SUCCESS;
}

// Entry point
int main(int argc, char* argv[])
{
    // TODO(bbowman)
    // Dispatch PacBio::CLI::Run on APPNAME for different commandline interfaces
    if (argc == 4) {
        return SimulateGenomicReads(argv[1], argv[2], argv[3]);
    } else {
        std::cerr << "genomic_sim takes exactly three arguments: <reference fasta> <input bam> "
                     "<output bam>\n";
        exit(EXIT_FAILURE);
    }
}
