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

#include <array>
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>

#include <pacbio/data/ArrayRead.h>
#include <pacbio/data/MSA.h>
#include <pacbio/io/BamParser.h>
#include <pacbio/io/Utility.h>
#include <pacbio/juliet/AminoAcidCaller.h>
#include <pacbio/juliet/JulietSettings.h>
#include <pacbio/juliet/ResistanceCaller.h>
#include <pacbio/statistics/Fisher.h>
#include <pacbio/statistics/Tests.h>
#include <pbbam/BamReader.h>
#include <pbbam/BamRecord.h>
#include <pbcopper/json/JSON.h>

#include <pacbio/juliet/JulietWorkflow.h>

namespace PacBio {
namespace Juliet {

std::ostream& JulietWorkflow::LogCI(const std::string& prefix)
{
    std::cout << std::setw(20) << std::left << prefix << ": ";
    return std::cout;
}

void JulietWorkflow::Run(const JulietSettings& settings)
{
    auto globalOutputPrefix = settings.OutputPrefix;
    globalOutputPrefix += globalOutputPrefix.empty() ? "" : "/";

    if (settings.Mode == AnalysisMode::BASE) {
        std::unordered_map<std::string, JSON::Json> jsonResults;
        for (const auto& inputFile : settings.InputFiles) {
            const auto outputPrefix = globalOutputPrefix + IO::FilePrefix(inputFile);

            // Convert BamRecords to unrolled ArrayReads
            std::vector<Data::ArrayRead> reads;
            reads = IO::ParseBam(inputFile, settings.RegionStart, settings.RegionEnd);

            Data::MSA msa(reads);

            // Compute fisher's exact test for each position
            for (auto& column : msa) {
                column.AddFisherResult(Statistics::Tests::FisherCCS(column));
                column.AddFisherResult(Statistics::Tests::FisherCCS(column, column.insertions));
            }

            if (settings.SaveMSA) {
                // Store msa + p-values
                std::ofstream msaStream(outputPrefix + ".msa");
                msaStream << "pos A Fa C Fc G Fg T Ft N Fn" << std::endl;
                int pos = msa.beginPos;
                for (auto& column : msa)
                    msaStream << ++pos << " " << column << std::endl;
                msaStream.close();
            }
            ResistanceCaller resiCaller(msa);

            const auto json = resiCaller.JSON();
            jsonResults.insert({IO::FilePrefix(inputFile), json});
            std::ofstream jsonStream(outputPrefix + ".json");
            jsonStream << json.dump(2) << std::endl;

            std::ofstream htmlStream(outputPrefix + ".html");
            ResistanceCaller::HTML(htmlStream, json, settings.DRMOnly, settings.Details);
        }
    } else if (settings.Mode == AnalysisMode::AMINO) {
        for (const auto& inputFile : settings.InputFiles) {
            const auto outputPrefix = globalOutputPrefix + IO::FilePrefix(inputFile);

            ErrorEstimates error;
            if (settings.SubstitutionRate != 0.0 && settings.DeletionRate != 0.0) {
                error = ErrorEstimates(settings.SubstitutionRate, settings.DeletionRate);
            } else {
                std::string chemistry;
                const BAM::BamReader bamReader(inputFile);
                const auto readGroups = bamReader.Header().ReadGroups();
                for (const auto& rg : readGroups) {
                    if (chemistry.empty())
                        chemistry = rg.SequencingChemistry();
                    else if (chemistry != rg.SequencingChemistry())
                        throw std::runtime_error("Mixed chemistries are not allowed");
                }
                error = ErrorEstimates(chemistry);
            }

            // Convert BamRecords to unrolled ArrayReads
            std::vector<Data::ArrayRead> reads;
            reads = IO::ParseBam(inputFile, settings.RegionStart, settings.RegionEnd);

            // Call variants
            AminoAcidCaller aac(reads, error, settings.TargetConfigUser);
            const auto json = aac.JSON();

            std::ofstream jsonStream(outputPrefix + ".json");
            jsonStream << json.dump(2) << std::endl;

            std::ofstream htmlStream(outputPrefix + ".html");
            AminoAcidCaller::HTML(htmlStream, json, settings.DRMOnly, settings.Details);

            // Store msa + p-values
            if (settings.SaveMSA) {
                std::ofstream msaStream(outputPrefix + ".msa");
                msaStream << "pos A C G T N" << std::endl;
                int pos = aac.msa_->beginPos;
                for (auto& column : *aac.msa_) {
                    msaStream << ++pos;
                    const std::array<int, 5>& counts = column;
                    for (const auto& c : counts)
                        msaStream << " " << c;
                    msaStream << std::endl;
                }
                msaStream.close();
            }
        }
    } else if (settings.Mode == AnalysisMode::PHASING) {
        for (const auto& inputFile : settings.InputFiles) {
            const auto outputPrefix = globalOutputPrefix + IO::FilePrefix(inputFile);

            // Convert BamRecords to unrolled ArrayReads
            std::vector<Data::ArrayRead> reads;
            reads = IO::ParseBam(inputFile, settings.RegionStart, settings.RegionEnd);

            Data::MSA msa(reads);

            // Compute fisher's exact test for each position
            for (auto& column : msa) {
                column.AddFisherResult(Statistics::Tests::FisherCCS(column));
                column.AddFisherResult(Statistics::Tests::FisherCCS(column, column.insertions));
            }

            Data::MSA msaWithPrior(reads, msa);
        }
    } else if (settings.Mode == AnalysisMode::ERROR) {
        for (const auto& inputFile : settings.InputFiles) {
            std::vector<Data::ArrayRead> reads;
            reads = IO::ParseBam(inputFile, settings.RegionStart, settings.RegionEnd);
            Data::MSA msa(reads);
            double sub = 0;
            double del = 0;
            int columnCount = 0;
            for (const auto& column : msa) {
                if (column.Coverage() > 100) {
                    del += column.Frequency(4);
                    sub += 1.0 - column.Frequency(4) - column.Frequency(column.MaxElement());
                    ++columnCount;
                }
            }
            std::cout << inputFile << std::endl;
            std::cout << "sub: " << (sub / columnCount) << std::endl;
            std::cout << "del: " << (del / columnCount) << std::endl;
        }
    }
}
}
}  // ::PacBio::Juliet