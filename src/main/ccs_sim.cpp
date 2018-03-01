// Author: David Seifert

#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <pbbam/BamReader.h>
#include <pbbam/BamRecord.h>
#include <pbbam/BamRecordBuilder.h>
#include <pbbam/BamWriter.h>
#include <pbbam/Cigar.h>
#include <pbbam/ReadGroupInfo.h>
#include <pbbam/TagCollection.h>

#include <pacbio/UnanimityVersion.h>
#include <pacbio/consensus/ModelConfig.h>
#include "../ModelFactory.h"
#include "../Simulator.h"

using namespace PacBio::BAM;
using namespace PacBio::Consensus;
using namespace PacBio::Data;

// these strings are part of the BAM header, they CANNOT contain newlines
const static std::string DESCRIPTION("Simulate (sub)reads from templates.");
const static std::string APPNAME("ccs_sim");

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

static int SimulateReads(const std::string& inputFilename, const std::string& outputFilename,
                         unsigned int seed = 42)
{
    BamReader reader{inputFilename};
    BamRecord inputRecord;

    BamHeader newHeader{PrepareHeader("ccs_sim", reader.Header().ReadGroups())};
    assert(newHeader.ReadGroups().size() == 1);
    const std::string newRg{newHeader.ReadGroups().front().Id()};

    std::vector<BamRecord> newRecords;
    std::default_random_engine rng{seed};

    int32_t zmw = 0;
    while (reader.GetNext(inputRecord)) {
        ++zmw;

        // 1. add old CCS as reference
        const std::string ccsName{inputRecord.FullName()};
        const std::string ccsSeq{inputRecord.Sequence(Orientation::GENOMIC)};
        const ReadGroupInfo ccsRg{inputRecord.ReadGroup()};

        newHeader.AddSequence(SequenceInfo{ccsName, std::to_string(ccsSeq.length())});
        const int32_t ccsId = newHeader.SequenceId(ccsName);

        // 2. simulate new read
        std::unique_ptr<ModelConfig> currentModel =
            ModelFactory::Create(ccsRg.SequencingChemistry(), {0, 0, 0, 0});
        std::pair<Read, std::vector<MoveType>> rawRead =
            currentModel->SimulateRead(&rng, ccsSeq, "");

        // 3. prepare new subread
        BamRecord newRecord{newHeader};

        const Cigar newCigar{ConvertStatePathToCigar(rawRead.second, ccsSeq, rawRead.first.Seq)};

        newRecord.ReadGroup(newRg)
            .IPD(Frames::Decode(rawRead.first.IPD), FrameEncodingType::LOSSY)
            .NumPasses(1)
            .PulseWidth(Frames::Decode(rawRead.first.PulseWidth), FrameEncodingType::LOSSY)
            .QueryStart(0)
            .QueryEnd(rawRead.first.Seq.length())
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
            .ReferenceId(ccsId)
            .SetMapped(true)
            .SetSequenceAndQualities(rawRead.first.Seq);

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
    // TODO(dseifert)
    // Dispatch PacBio::CLI::Run on APPNAME for different commandline interfaces
    if (argc == 3) {
        return SimulateReads(argv[1], argv[2]);
    } else {
        std::cerr << "ccs_sim takes exactly two arguments: <input bam> <output bam>\n";
        exit(EXIT_FAILURE);
    }
}
