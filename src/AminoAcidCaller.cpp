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
#include <functional>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <unordered_map>
#include <vector>

#include <pacbio/data/MSA.h>
#include <pacbio/juliet/AminoAcidCaller.h>
#include <pacbio/statistics/Fisher.h>
#include <pbcopper/json/JSON.h>

namespace PacBio {
namespace Juliet {
AminoAcidCaller::AminoAcidCaller(const std::vector<Data::ArrayRead>& reads,
                                 const ErrorModel& errorModel, const TargetConfig& targetConfig)
    : errorModel_(errorModel), targetConfig_(targetConfig)
{
    for (const auto& r : reads) {
        beginPos_ = std::min(beginPos_, r.ReferenceStart());
        endPos_ = std::max(endPos_, r.ReferenceEnd());
    }
    msa_ = std::make_unique<Data::MSA>(reads);

    GenerateMSA(reads);

    beginPos_ += 1;
    endPos_ += 1;

    CallVariants(reads);
}

void AminoAcidCaller::GenerateMSA(const std::vector<Data::ArrayRead>& reads)
{
    matrix_.reserve(reads.size());
    for (const auto& r : reads) {
        int pos = r.ReferenceStart() - beginPos_;
        assert(pos >= 0);
        std::vector<char> row(endPos_ - beginPos_, ' ');
        for (const auto& b : r.Bases) {
            switch (b.Cigar) {
                case 'X':
                case '=':
                    row[pos++] = b.Nucleotide;
                    break;
                case 'D':
                    row[pos++] = '-';
                    break;
                case 'I':
                case 'P':
                    break;
                default:
                    throw std::runtime_error("Unexpected cigar " + std::to_string(b.Cigar));
            }
        }
        matrix_.emplace_back(row);
    }
}

int AminoAcidCaller::CountNumberOfTests(const std::vector<TargetGene>& genes) const
{
    int numberOfTests = 0;
    for (const auto& gene : genes) {
        for (int i = gene.begin; i < gene.end - 2; ++i) {
            // Relative to gene begin
            const int ri = i - gene.begin;
            // Only work on beginnings of a codon
            if (ri % 3 != 0) continue;
            // Relative to window begin
            const int bi = i - beginPos_;

            std::unordered_map<std::string, int> codons;
            int coverage = 0;
            for (const auto& row : matrix_) {
                // Read does not cover codon
                if (bi + 2 >= static_cast<int>(row.size()) || bi < 0) continue;
                if (row.at(bi + 0) == ' ' || row.at(bi + 1) == ' ' || row.at(bi + 2) == ' ')
                    continue;
                ++coverage;

                // Read has a deletion
                if (row.at(bi + 0) == '-' || row.at(bi + 1) == '-' || row.at(bi + 2) == '-')
                    continue;

                std::string codon = std::string() + row.at(bi) + row.at(bi + 1) + row.at(bi + 2);

                // Codon is bogus
                if (codonToAmino_.find(codon) == codonToAmino_.cend()) continue;

                codons[codon]++;
            }
            numberOfTests += codons.size();
        }
    }
    return numberOfTests;
}

std::string AminoAcidCaller::FindDRMs(const std::string& geneName,
                                      const std::vector<TargetGene>& genes,
                                      const int position) const
{
    std::string drmSummary;
    for (const auto& gene : genes) {
        if (geneName == gene.name) {
            for (const auto& drms : gene.drms) {
                if (std::find(drms.positions.cbegin(), drms.positions.cend(), position) !=
                    drms.positions.cend()) {
                    if (!drmSummary.empty()) drmSummary += " + ";
                    drmSummary += drms.name;
                }
            }
            break;
        }
    }
    return drmSummary;
};
void AminoAcidCaller::CallVariants(const std::vector<Data::ArrayRead>& reads)
{
    auto genes = targetConfig_.targetGenes;
    bool hasReference = !targetConfig_.referenceSequence.empty();
    // If no user config has been provided, use complete input region
    if (genes.empty()) {
        TargetGene tg(beginPos_, endPos_, "unknown", {});
        genes.emplace_back(tg);
    }

    const ErrorEstimates error(errorModel_);
    VariantGene curVariantGene;
    std::string geneName;
    int geneOffset = 0;

    const auto SetNewGene = [this, &geneName, &curVariantGene, &geneOffset](
        const int begin, const std::string& name) {
        geneName = name;
        if (!curVariantGene.relPositionToVariant.empty())
            variantGenes_.push_back(std::move(curVariantGene));
        curVariantGene = VariantGene();
        curVariantGene.geneName = name;
        geneOffset = begin;
    };

    const auto CodonProbability = [&error](const std::string& a, const std::string& b) {
        double p = 1;
        for (int i = 0; i < 3; ++i) {
            if (a[i] == '-' || b[i] == '-')
                p *= error.deletion;
            else if (a[i] != b[i])
                p *= error.substitution;
            else
                p *= error.match;
        }
        return p;
    };

    const int numberOfTests = CountNumberOfTests(genes);

#ifdef JULIET_INHOUSE_PERFORMANCE
    double truePositives = 0;
    double falsePositives = 0;
    double falseNegative = 0;
    double trueNegative = 0;
    auto MeasurePerformance = [&truePositives, &falsePositives, &falseNegative, &trueNegative, this,
                               &geneName](const auto& codon_counts, const auto& codonPos,
                                          const auto& i, const auto& p) {
        const auto curCodon = codonToAmino_.at(codon_counts.first);
        bool predictor = (i == 3191 && curCodon == 'Y' && "TAC" == codon_counts.first) ||
                         (i == 2741 && curCodon == 'R' && "AGA" == codon_counts.first) ||
                         (i == 2669 && curCodon == 'L' && "TTG" == codon_counts.first);
        bool ignored =
            (geneName == "Protease" && codonPos == 3 && curCodon == 'I') ||
            (geneName == "Protease" && codonPos == 37 && curCodon == 'N') ||
            (geneName == "Reverse Transcriptase" && codonPos == 102 && curCodon == 'Q') ||
            (geneName == "Reverse Transcriptase" && codonPos == 122 && curCodon == 'K') ||
            (geneName == "Reverse Transcriptase" && codonPos == 162 && curCodon == 'C') ||
            (geneName == "Reverse Transcriptase" && codonPos == 214 && curCodon == 'F') ||
            (geneName == "Reverse Transcriptase" && codonPos == 272 && curCodon == 'A') ||
            (geneName == "Reverse Transcriptase" && codonPos == 293 && curCodon == 'V') ||
            (geneName == "Reverse Transcriptase" && codonPos == 358 && curCodon == 'K') ||
            (geneName == "Reverse Transcriptase" && codonPos == 376 && curCodon == 'A') ||
            (geneName == "Reverse Transcriptase" && codonPos == 400 && curCodon == 'A') ||
            (geneName == "Reverse Transcriptase" && codonPos == 435 && curCodon == 'I') ||
            (geneName == "RNase" && codonPos == 20 && curCodon == 'D') ||
            (geneName == "RNase" && codonPos == 28 && curCodon == 'P') ||
            (geneName == "RNase" && codonPos == 43 && curCodon == 'H') ||
            (geneName == "RNase" && codonPos == 72 && curCodon == 'K') ||
            (geneName == "RNase" && codonPos == 79 && curCodon == 'S') ||
            (geneName == "Integrase" && codonPos == 10 && curCodon == 'E') ||
            (geneName == "Integrase" && codonPos == 113 && curCodon == 'V') ||
            (geneName == "Integrase" && codonPos == 123 && curCodon == 'S') ||
            (geneName == "Integrase" && codonPos == 124 && curCodon == 'T') ||
            (geneName == "Integrase" && codonPos == 127 && curCodon == 'K') ||
            (geneName == "Integrase" && codonPos == 151 && curCodon == 'I') ||
            (geneName == "Integrase" && codonPos == 232 && curCodon == 'D') ||
            (geneName == "Integrase" && codonPos == 234 && curCodon == 'V');

        if (!ignored) {
            if (p < alpha) {
                if (predictor)
                    ++truePositives;
                else
                    ++falsePositives;
            } else {
                if (predictor)
                    ++falseNegative;
                else
                    ++trueNegative;
            }
        }

        return !ignored;
    };
#endif

    for (const auto& gene : genes) {
        SetNewGene(gene.begin, gene.name);
        for (int i = gene.begin; i < gene.end - 2; ++i) {
            // Absolute reference position
            const int ai = i - 1;
            // Relative to gene begin
            const int ri = i - geneOffset;
            // Only work on beginnings of a codon
            if (ri % 3 != 0) continue;
            // Relative to window begin
            const int bi = i - beginPos_;

            const int codonPos = 1 + (ri / 3);
            auto& curVariantPosition = curVariantGene.relPositionToVariant[codonPos];

            std::map<std::string, int> codons;
            int coverage = 0;
            for (const auto& row : matrix_) {
                // Read does not cover codon
                if (bi + 2 > static_cast<int>(row.size()) || bi < 0) continue;
                if (row.at(bi + 0) == ' ' || row.at(bi + 1) == ' ' || row.at(bi + 2) == ' ')
                    continue;
                ++coverage;

                // Read has a deletion
                if (row.at(bi + 0) == '-' || row.at(bi + 1) == '-' || row.at(bi + 2) == '-')
                    continue;

                const auto codon = std::string() + row.at(bi) + row.at(bi + 1) + row.at(bi + 2);

                // Codon is bogus
                if (codonToAmino_.find(codon) == codonToAmino_.cend()) continue;

                codons[codon]++;
            }

            if (hasReference) {
                curVariantPosition.refCodon = targetConfig_.referenceSequence.substr(ai, 3);
                if (codonToAmino_.find(curVariantPosition.refCodon) == codonToAmino_.cend()) {
                    continue;
                }
                curVariantPosition.refAminoAcid = codonToAmino_.at(curVariantPosition.refCodon);
            } else {
                int max = -1;
                std::string argmax;
                for (const auto& codon_counts : codons) {
                    if (codon_counts.second > max) {
                        max = codon_counts.second;
                        argmax = codon_counts.first;
                    }
                }
                curVariantPosition.refCodon = argmax;
                if (codonToAmino_.find(curVariantPosition.refCodon) == codonToAmino_.cend()) {
                    continue;
                }
                curVariantPosition.refAminoAcid = codonToAmino_.at(curVariantPosition.refCodon);
            }

            for (const auto& codon_counts : codons) {
                if (codonToAmino_.at(codon_counts.first) == curVariantPosition.refAminoAcid)
                    continue;
                double p = (Statistics::Fisher::fisher_exact_tiss(
                                codon_counts.second, coverage,
                                coverage * CodonProbability(curVariantPosition.refCodon,
                                                            codon_counts.first),
                                coverage) *
                            numberOfTests);

                if (p > 1) p = 1;

#ifdef JULIET_INHOUSE_PERFORMANCE
                const bool variableSite = MeasurePerformance(codon_counts, codonPos, ai, p);

                if (variableSite && p < alpha) {
#else
                if (p < alpha) {
#endif
                    VariantGene::VariantPosition::VariantCodon curVariantCodon;
                    curVariantCodon.codon = codon_counts.first;
                    curVariantCodon.frequency = codon_counts.second / static_cast<double>(coverage);
                    curVariantCodon.pValue = p;
                    curVariantCodon.knownDRM = FindDRMs(geneName, genes, codonPos);

                    curVariantPosition.aminoAcidToCodons[codonToAmino_.at(codon_counts.first)]
                        .push_back(curVariantCodon);
                }
            }
            if (!curVariantPosition.aminoAcidToCodons.empty()) {
                curVariantPosition.coverage = coverage;
                for (int j = -3; j < 6; ++j) {
                    if (i + j >= beginPos_ && i + j < endPos_) {
                        int abs = ai + j;
                        JSON::Json msaCounts;
                        msaCounts["rel_pos"] = j;
                        msaCounts["abs_pos"] = abs;
                        msaCounts["A"] = (*msa_)[abs][0];
                        msaCounts["C"] = (*msa_)[abs][1];
                        msaCounts["G"] = (*msa_)[abs][2];
                        msaCounts["T"] = (*msa_)[abs][3];
                        msaCounts["-"] = (*msa_)[abs][4];
                        if (hasReference)
                            msaCounts["wt"] =
                                std::string(1, targetConfig_.referenceSequence.at(abs));
                        else
                            msaCounts["wt"] =
                                std::string(1, Data::TagToNucleotide((*msa_)[abs].MaxElement()));
                        curVariantPosition.msa.push_back(msaCounts);
                    }
                }
            }
        }
    }
#ifdef JULIET_INHOUSE_PERFORMANCE
    std::cerr << (truePositives / 3.0);
    std::cerr << " " << (falsePositives / (numberOfTests - 3));
    std::cerr << " " << numberOfTests;
    std::cerr << " " << ((truePositives + trueNegative) /
                         (truePositives + falsePositives + falseNegative + trueNegative))
              << std::endl;
#endif
    if (!curVariantGene.relPositionToVariant.empty())
        variantGenes_.push_back(std::move(curVariantGene));
}

void AminoAcidCaller::HTML(std::ostream& out, const JSON::Json& j, bool onlyKnownDRMs, bool details)
{
#if 1
    auto strip = [](const std::string& input) -> std::string {
        std::string s = input;
        s.erase(std::remove(s.begin(), s.end(), '\"'), s.end());
        return s;
    };
    out << "<html>" << std::endl
        << "<head>" << std::endl
        << R"(
            <script src="http://ajax.googleapis.com/ajax/libs/jquery/1.11.1/jquery.min.js"></script>
            <script type="text/javascript">
            $(document).ready(function() {
                $(".var").bind( "click", function( event ) {
                    $(this).next().slideToggle(0);
            });
            });
            </script>)"
        << std::endl
        << "<style>" << std::endl
        << R"(
            body { font-family: helvetica-light }
            table { border-collapse: collapse; margin-bottom: 20px; }
            tr:nth-child(1) { background-color: #3d3d3d; color: white; }
            tr:nth-child(3) th { padding: 5px 5px 5px 5px; text-align: center; border-bottom: 1px solid #2d2d2d; }
            tr:nth-child(2) th:nth-child(2) { border-left: 1px dashed black; }
            tr:nth-child(3) th:nth-child(3) { border-right: 1px dashed black; }
            td { padding: 15px 5px 15px 5px; text-align: center; border-bottom: 1px solid white; }
            table td:nth-child(1) { background-color:#ddd; border-right: 1px solid #eee; }
            table td:nth-child(2) { background-color:#eee; border-right: 1px solid #ddd; }
            table td:nth-child(3) { background-color:#fff; border-right: 1px solid #ddd; font-weight: bold;}
            table td:nth-child(4) { background-color:#eee; border-right: 1px dashed #ccc;  }
            table td:nth-child(5) { background-color: #ddd; border-right: 1px dashed #bbb; }
            table td:nth-child(6) { background-color: #ccc; border-right: 1px dashed #aaa; }
            table td:nth-child(7) { background-color: #bbb;}
            table td:nth-child(8) { background-color: #aaa; color: #fff600}
            tr:not(.msa):hover td { background-color: white; }
            tr:not(.msa):hover td:nth-child(8) { color: purple; }
            .msa table tr:hover td { background-color: gray; color:white; }
            .top table { background-color:white; border:0; }
            .top table td { background-color:white; border:0; border-bottom: 1px solid gray; font-weight: normal}
            .top table tr { border:0; }
            .top table th { border:0; }
            .msa { display:none; }
            )"
        << std::endl
        << "</style>" << std::endl
        << "</head>" << std::endl
        << R"(<body>
            <details style="margin-bottom: 20px">
            <summary>Legend</summary>
            <p>Every table represents a gene in the Pol polyprotein.<br/>
            Each row stands for a mutated amino acid. Positions are relative to the current gene.<br/>
            Positions with no or synonymous mutation are not being shown.<br/>
            The used reference is HXB2 and all coordinates are in reference space.<br/>
            The mutated nucleotide is highlighted in the codon.<br/>
            Frequency is per codon.<br/>
            Coverage includes deletions.<br/>
            Known drug-resistance mutations positions are annotated in the last column,<br/>
            whereas 'S' stands for surveillance. Annotations from the <a href="https://hivdb.stanford.edu" target="_new">Stanford DB</a>.<br/>
            <br/>
            Clicking on a row unfolds the counts of the multiple sequence alignment of the<br/>
            codon position and up to +-3 surrounding positions.<br/>
            Red colored are nucleotides of the codon and in bold the wild type.<br/>
            <br/>
            Deletions and insertions are being ignored in this version.<br/>
            <br/>
            This software is for research only and has not been clinically validated!</p>
            </details>)"
        << std::endl;

    if (j.find("genes") == j.cend() || j["genes"].is_null()) return;
    for (const auto& gene : j["genes"]) {
        out << "<table class=\"top\">" << std::endl
            << R"(
                <col width="40px"/>
                <col width="40px"/>
                <col width="40px"/>
                <col width="40px"/>
                <col width="40px"/>
                <col width="60px"/>
                <col width="60px"/>
                <col width="180px"/>
                <tr>
                <th colspan="9">)"
            << strip(gene["name"]) << R"(</th>
                </tr>
                <tr>
                <th colspan="3">HXB2</th>
                <th colspan="5">Sample</th>
                </tr>
                <tr>
                <th>Codon</th>
                <th>AA</th>
                <th>Pos</th>
                <th>AA</th>
                <th colspan="1">Codon</th>
                <th colspan="1">Frequency</th>
                <th colspan="1">Coverage</th>
                <th colspan="1">DRM</th>
                </tr>)"
            << std::endl;

        for (auto& variantPosition : gene["variant_positions"]) {
            std::stringstream line;
            const std::string refCodon = strip(variantPosition["ref_codon"]);
            line << "<tr class=\"var\">\n"
                 << "<td>" << strip(variantPosition["ref_codon"])[0] << " "
                 << strip(variantPosition["ref_codon"])[1] << " "
                 << strip(variantPosition["ref_codon"])[2] << "</td>\n"
                 << "<td>" << strip(variantPosition["ref_amino_acid"]) << "</td>\n"
                 << "<td>" << variantPosition["ref_position"] << "</td>";
            std::string prefix = line.str();
            line.str("");
            bool first = true;
            for (auto& variant_amino_acid : variantPosition["variant_amino_acids"]) {
                for (auto& variant_codons : variant_amino_acid["variant_codons"]) {
                    bool mutated[]{
                        strip(variantPosition["ref_codon"])[0] != strip(variant_codons["codon"])[0],
                        strip(variantPosition["ref_codon"])[1] != strip(variant_codons["codon"])[1],
                        strip(variantPosition["ref_codon"])[2] !=
                            strip(variant_codons["codon"])[2]};
                    line << "<td>" << strip(variant_amino_acid["amino_acid"]) << "</td>";
                    line << "<td>";
                    for (int j = 0; j < 3; ++j) {
                        if (mutated[j]) line << "<b style=\"color:#ff5e5e; font-weight:normal\">";
                        line << strip(variant_codons["codon"])[j] << " ";
                        if (mutated[j]) line << "</b>";
                    }

                    double fOrig = variant_codons["frequency"];
                    double fTmp;
                    int exp = 0;
                    do {
                        fTmp = fOrig * std::pow(10, ++exp);
                    } while (static_cast<int>(fTmp) < 10);
                    fOrig = static_cast<int>(fOrig * std::pow(10, exp));
                    fOrig /= std::pow(10, exp);
                    line << "<td>" << fOrig << "</td>";
                    if (first) {
                        out << prefix << line.str();
                        out << "<td>" << variantPosition["coverage"] << "</td>";
                        first = false;
                    } else {
                        out << "<tr class=\"var\"><td></td><td></td><td></td>" << line.str()
                            << "<td></td>";
                    }
                    out << "<td>" << strip(variant_codons["known_drm"]) << "</td>";
                    out << "</tr>" << std::endl;
                    line.str("");

                    out << R"(
                        <tr class="msa">
                        <td colspan=3 style="background-color: white"></td>
                        <td colspan=14 style="padding:0; margin:0">
                        <table style="padding:0; margin:0">
                        <col width="80px" />
                        <col width="80px" />
                        <col width="80px" />
                        <col width="80px" />
                        <col width="80px" />
                        <col width="80px" />
                        <tr style="padding:0">
                        <th style="padding:2px 0 0px 0">Pos</th>
                        <th style="padding:2px 0 0px 0">A</th>
                        <th style="padding:2px 0 0px 0">C</th>
                        <th style="padding:2px 0 0px 0">G</th>
                        <th style="padding:2px 0 0px 0">T</th>
                        <th style="padding:2px 0 0px 0">-</th>
                        </tr>
                        )";

                    for (auto& column : variantPosition["msa"]) {
                        int relPos = column["rel_pos"];
                        out << "<tr><td>" << relPos << "</td>" << std::endl;
                        for (int j = 0; j < 5; ++j) {
                            out << "<td style=\"";
                            if (relPos >= 0 && relPos < 3) {
                                if (j ==
                                    Data::NucleotideToTag(strip(variant_codons["codon"])[relPos]))
                                    out << "color:red;";
                            }
                            if (j == Data::NucleotideToTag(strip(column["wt"])[0]))
                                out << "font-weight:bold;";
                            out << "\">" << column[std::string(1, Data::TagToNucleotide(j))]
                                << "</td>" << std::endl;
                        }
                        out << "</tr>" << std::endl;
                    }
                    out << "</table></tr>" << std::endl;
                }
            }
        }
    }
    out << "</table>" << std::endl << "</body>" << std::endl << "</html>" << std::endl;
#endif
}

JSON::Json AminoAcidCaller::JSON()
{
    using namespace JSON;
    Json root;
    std::vector<Json> genes;
    for (const auto& v : variantGenes_) {
        Json j = v.ToJson();
        if (j.find("variant_positions") != j.cend()) genes.push_back(j);
    }
    root["genes"] = genes;

    return root;
}

const std::unordered_map<std::string, char> AminoAcidCaller::codonToAmino_ = {
    {"ATT", 'I'}, {"ATC", 'I'}, {"ATA", 'I'}, {"CTT", 'L'}, {"CTC", 'L'}, {"CTA", 'L'},
    {"CTG", 'L'}, {"TTA", 'L'}, {"TTG", 'L'}, {"GTT", 'V'}, {"GTC", 'V'}, {"GTA", 'V'},
    {"GTG", 'V'}, {"TTT", 'F'}, {"TTC", 'F'}, {"ATG", 'M'}, {"TGT", 'C'}, {"TGC", 'C'},
    {"GCT", 'A'}, {"GCC", 'A'}, {"GCA", 'A'}, {"GCG", 'A'}, {"GGT", 'G'}, {"GGC", 'G'},
    {"GGA", 'G'}, {"GGG", 'G'}, {"CCT", 'P'}, {"CCC", 'P'}, {"CCA", 'P'}, {"CCG", 'P'},
    {"ACT", 'T'}, {"ACC", 'T'}, {"ACA", 'T'}, {"ACG", 'T'}, {"TCT", 'S'}, {"TCC", 'S'},
    {"TCA", 'S'}, {"TCG", 'S'}, {"AGT", 'S'}, {"AGC", 'S'}, {"TAT", 'Y'}, {"TAC", 'Y'},
    {"TGG", 'W'}, {"CAA", 'Q'}, {"CAG", 'Q'}, {"AAT", 'N'}, {"AAC", 'N'}, {"CAT", 'H'},
    {"CAC", 'H'}, {"GAA", 'E'}, {"GAG", 'E'}, {"GAT", 'D'}, {"GAC", 'D'}, {"AAA", 'K'},
    {"AAG", 'K'}, {"CGT", 'R'}, {"CGC", 'R'}, {"CGA", 'R'}, {"CGG", 'R'}, {"AGA", 'R'},
    {"AGG", 'R'}, {"TAA", 'X'}, {"TAG", 'X'}, {"TGA", 'X'}};
}
}  // ::PacBio::Juliet
