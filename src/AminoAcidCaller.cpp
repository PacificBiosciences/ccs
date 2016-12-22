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
                                 const ErrorModel& errorModel)
    : errorModel_(errorModel)
{
    for (const auto& r : reads) {
        beginPos_ = std::min(beginPos_, r.ReferenceStart());
        endPos_ = std::max(endPos_, r.ReferenceEnd());
    }
    msa_ = std::make_unique<Data::MSA>(reads);

    GenerateMSA(reads);
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

struct Gene
{
    int begin;
    int end;
    std::string name;
};

void AminoAcidCaller::CallVariants(const std::vector<Data::ArrayRead>& reads)
{
    std::vector<Gene> genes{
        {2253, 2550, "Protease"},
        {2550, 3870, "Reverse Transcriptase"},
        {3870, 4230, "RNase"},
        {4230, 5096, "Integrase"}
    };

    ErrorEstimates error(errorModel_);
    VariantGene curVariantGene;
    std::string geneName;
    int geneOffset = 0;
    auto SetNewGene = [this, &geneName, &curVariantGene, &geneOffset](const int begin,
                                                                      const std::string& name) {
        geneName = name;
        if (!curVariantGene.relPositionToVariant.empty())
            variantGenes_.push_back(std::move(curVariantGene));
        curVariantGene = VariantGene();
        curVariantGene.geneName = name;
        geneOffset = begin;
    };
    auto CodonProbability = [&error](const std::string& a, const std::string& b) {
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
    auto FindDRMs = [this, &geneName](int position) {
        std::string drm;
        if (geneName == "Protease") {
            if (std::find(pi.cbegin(), pi.cend(), position) != pi.cend())
                drm = "PI";
            else if (std::find(piSurveillance.cbegin(), piSurveillance.cend(), position) !=
                     piSurveillance.cend())
                drm = "PI S";
        } else if (geneName == "Reverse Transcriptase") {
            if (std::find(nnrti.cbegin(), nnrti.cend(), position) != nnrti.cend())
                drm = "NNRTI";
            else if (std::find(nnrtiSurveillance.cbegin(), nnrtiSurveillance.cend(), position) !=
                     nnrtiSurveillance.cend())
                drm = "NNRTI S";

            if (std::find(nrti.cbegin(), nrti.cend(), position) != nrti.cend()) {
                if (!drm.empty()) drm += " + ";
                drm += "NRTI";
            } else if (std::find(nrtiSurveillance.cbegin(), nrtiSurveillance.cend(), position) !=
                       nrtiSurveillance.cend()) {
                if (!drm.empty()) drm += " + ";
                drm += "NRTI S";
            }
        } else if (geneName == "Integrase") {
            if (std::find(ini.cbegin(), ini.cend(), position) != ini.cend()) drm = "INSTI";
        }
        return drm;
    };
    int numberOfTests = 0;
    for (int i = beginPos_; i < endPos_ - 2; ++i) {
        // Corrected index for 1-based reference
        const int ci = i + 1;
        bool missing = true;
        for (const auto& gene : genes)
        {
            if (ci >= gene.begin && ci < gene.end)
            {
                geneOffset = gene.begin;
                break;
            }
        }
        if (ci >= genes.back().end) geneOffset = genes.back().end;

        // Relative to gene begin
        int ri = ci - geneOffset;
        // Only work on beginnings of a codon
        if (ri % 3 != 0) continue;
        // Relative to window begin
        int bi = i - beginPos_;

        int codonPos = (1 + ri / 3);
        auto& curVariantPosition = curVariantGene.relPositionToVariant[codonPos];

        curVariantPosition.refCodon = ref_.substr(ci - 1, 3);
        if (codonToAmino_.find(curVariantPosition.refCodon) == codonToAmino_.cend()) {
            continue;
        }

        std::map<std::string, int> codons;
        int coverage = 0;
        for (const auto& row : matrix_) {
            // Read does not cover codon
            if (row.at(bi + 0) == ' ' || row.at(bi + 1) == ' ' || row.at(bi + 2) == ' ') continue;
            ++coverage;

            // Read has a deletion
            if (row.at(bi + 0) == '-' || row.at(bi + 1) == '-' || row.at(bi + 2) == '-') continue;

            std::string codon = std::string() + row.at(bi) + row.at(bi + 1) + row.at(bi + 2);

            // Codon is bogus
            if (codonToAmino_.find(codon) == codonToAmino_.cend()) continue;

            codons[codon]++;
        }
        numberOfTests += codons.size();
    }

    geneOffset = 0;

#ifdef PERFORMANCE
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
            if (p < 0.01) {
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

    for (int i = beginPos_; i < endPos_ - 2; ++i) {
        // Corrected index for 1-based reference
        const int ci = i + 1;

        bool missing = true;
        for (const auto& gene : genes)
        {
            if (ci >= gene.begin && ci < gene.end && geneName != gene.name)
            {
                SetNewGene(gene.begin, gene.name);
                break;
            }
        }
        if (ci >= genes.back().end) {
            SetNewGene(genes.back().end, genes.back().name);
        }

        // Relative to gene begin
        int ri = ci - geneOffset;
        // Only work on beginnings of a codon
        if (ri % 3 != 0) continue;
        // Relative to window begin
        int bi = i - beginPos_;

        int codonPos = (1 + ri / 3);
        auto& curVariantPosition = curVariantGene.relPositionToVariant[codonPos];

        curVariantPosition.refCodon = ref_.substr(ci - 1, 3);
        if (codonToAmino_.find(curVariantPosition.refCodon) == codonToAmino_.cend()) {
            continue;
        }
        curVariantPosition.refAminoAcid = codonToAmino_.at(curVariantPosition.refCodon);

        std::map<std::string, int> codons;
        int coverage = 0;
        for (const auto& row : matrix_) {
            // Read does not cover codon
            if (row.at(bi + 0) == ' ' || row.at(bi + 1) == ' ' || row.at(bi + 2) == ' ') continue;
            ++coverage;

            // Read has a deletion
            if (row.at(bi + 0) == '-' || row.at(bi + 1) == '-' || row.at(bi + 2) == '-') continue;

            std::string codon = std::string() + row.at(bi) + row.at(bi + 1) + row.at(bi + 2);

            // Codon is bogus
            if (codonToAmino_.find(codon) == codonToAmino_.cend()) continue;

            codons[codon]++;
        }
        for (const auto& codon_counts : codons) {
            if (codonToAmino_.at(codon_counts.first) == curVariantPosition.refAminoAcid) continue;
            double p =
                (Statistics::Fisher::fisher_exact_tiss(
                     codon_counts.second, coverage,
                     coverage * CodonProbability(curVariantPosition.refCodon, codon_counts.first),
                     coverage) *
                 numberOfTests);

            if (p > 1) p = 1;

#ifdef PERFORMANCE
            bool variableSite = MeasurePerformance(codon_counts, codonPos, i, p);

            if (variableSite && p < 0.01) {
#else
            if (p < 0.01) {
#endif
                VariantGene::VariantPosition::VariantCodon curVariantCodon;
                curVariantCodon.codon = codon_counts.first;
                curVariantCodon.frequency = codon_counts.second / static_cast<double>(coverage);
                curVariantCodon.pValue = p;
                curVariantCodon.knownDRM = FindDRMs(codonPos);

                curVariantPosition.aminoAcidToCodons[codonToAmino_.at(codon_counts.first)]
                    .push_back(curVariantCodon);
            }
        }
        if (!curVariantPosition.aminoAcidToCodons.empty()) {
            curVariantPosition.coverage = coverage;
            for (int j = -3; j < 6; ++j) {
                if (i + j >= beginPos_ && i + j < endPos_) {
                    JSON::Json msaCounts;
                    msaCounts["rel_pos"] = j;
                    msaCounts["abs_pos"] = i + j;
                    msaCounts["A"] = (*msa_)[i + j][0];
                    msaCounts["C"] = (*msa_)[i + j][1];
                    msaCounts["G"] = (*msa_)[i + j][2];
                    msaCounts["T"] = (*msa_)[i + j][3];
                    msaCounts["-"] = (*msa_)[i + j][4];
                    msaCounts["wt"] = std::string(1, ref_[i + j]);
                    curVariantPosition.msa.push_back(msaCounts);
                }
            }
        }
    }
#ifdef PERFORMANCE
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
    out << "<html>" << std::endl;
    out << "<head>" << std::endl;
    out << R"(
<script src="http://ajax.googleapis.com/ajax/libs/jquery/1.11.1/jquery.min.js"></script>
<script type="text/javascript">
$(document).ready(function() {
    $(".var").bind( "click", function( event ) {
        $(this).next().slideToggle(0);
});
});
</script>)"
        << std::endl;
    out << "<style>" << std::endl;
    out << R"(
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
)" << std::endl;
    out << "</style>" << std::endl;
    out << "</head>" << std::endl;
    out << R"(<body>
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
        out << "<table class=\"top\">" << std::endl;
        out << R"(
<col width="40px"/>
<col width="40px"/>
<col width="40px"/>
<col width="40px"/>
<col width="40px"/>
<col width="60px"/>
<col width="60px"/>
<col width="180px"/>
<tr>
<th colspan="9">)";
        out << strip(gene["name"]);
        out << R"(</th>
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
</tr>)" << std::endl;

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
    out << "</table>" << std::endl;
    out << "</body>" << std::endl;
    out << "</html>" << std::endl;
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

const std::string AminoAcidCaller::ref_ =
    "TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCTGATTAGCAG"
    "AACTACACACCAGGGCCAGGGATCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTACCAGTTGAGCCAGAGAAGTTAGA"
    "AGAAGCCAACAAAGGAGAGAACACCAGCTTGTTACACCCTGTGAGCCTGCATGGAATGGATGACCCGGAGAGAGAAGTGTTAGAGT"
    "GGAGGTTTGACAGCCGCCTAGCATTTCATCACATGGCCCGAGAGCTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGC"
    "TACAAGGGACTTTCCGCTGGGGACTTTCCAGGGAGGCGTGGCCTGGGCGGGACTGGGGAGTGGCGAGCCCTCAGATCCTGCATATA"
    "AGCAGCTGCTTTTTGCCTGTACTGGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTG"
    "CTTAAGCCTCAATAAAGCTTGCCTTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTCAG"
    "ACCCTTTTAGTCAGTGTGGAAAATCTCTAGCAGTGGCGCCCGAACAGGGACCTGAAAGCGAAAGGGAAACCAGAGGAGCTCTCTCG"
    "ACGCAGGACTCGGCTTGCTGAAGCGCGCACGGCAAGAGGCGAGGGGCGGCGACTGGTGAGTACGCCAAAAATTTTGACTAGCGGAG"
    "GCTAGAAGGAGAGAGATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGG"
    "GGGAAAGAAAAAATATAAATTAAAACATATAGTATGGGCAAGCAGGGAGCTAGAACGATTCGCAGTTAATCCTGGCCTGTTAGAAA"
    "CATCAGAAGGCTGTAGACAAATACTGGGACAGCTACAACCATCCCTTCAGACAGGATCAGAAGAACTTAGATCATTATATAATACA"
    "GTAGCAACCCTCTATTGTGTGCATCAAAGGATAGAGATAAAAGACACCAAGGAAGCTTTAGACAAGATAGAGGAAGAGCAAAACAA"
    "AAGTAAGAAAAAAGCACAGCAAGCAGCAGCTGACACAGGACACAGCAATCAGGTCAGCCAAAATTACCCTATAGTGCAGAACATCC"
    "AGGGGCAAATGGTACATCAGGCCATATCACCTAGAACTTTAAATGCATGGGTAAAAGTAGTAGAAGAGAAGGCTTTCAGCCCAGAA"
    "GTGATACCCATGTTTTCAGCATTATCAGAAGGAGCCACCCCACAAGATTTAAACACCATGCTAAACACAGTGGGGGGACATCAAGC"
    "AGCCATGCAAATGTTAAAAGAGACCATCAATGAGGAAGCTGCAGAATGGGATAGAGTGCATCCAGTGCATGCAGGGCCTATTGCAC"
    "CAGGCCAGATGAGAGAACCAAGGGGAAGTGACATAGCAGGAACTACTAGTACCCTTCAGGAACAAATAGGATGGATGACAAATAAT"
    "CCACCTATCCCAGTAGGAGAAATTTATAAAAGATGGATAATCCTGGGATTAAATAAAATAGTAAGAATGTATAGCCCTACCAGCAT"
    "TCTGGACATAAGACAAGGACCAAAGGAACCCTTTAGAGACTATGTAGACCGGTTCTATAAAACTCTAAGAGCCGAGCAAGCTTCAC"
    "AGGAGGTAAAAAATTGGATGACAGAAACCTTGTTGGTCCAAAATGCGAACCCAGATTGTAAGACTATTTTAAAAGCATTGGGACCA"
    "GCGGCTACACTAGAAGAAATGATGACAGCATGTCAGGGAGTAGGAGGACCCGGCCATAAGGCAAGAGTTTTGGCTGAAGCAATGAG"
    "CCAAGTAACAAATTCAGCTACCATAATGATGCAGAGAGGCAATTTTAGGAACCAAAGAAAGATTGTTAAGTGTTTCAATTGTGGCA"
    "AAGAAGGGCACACAGCCAGAAATTGCAGGGCCCCTAGGAAAAAGGGCTGTTGGAAATGTGGAAAGGAAGGACACCAAATGAAAGAT"
    "TGTACTGAGAGACAGGCTAATTTTTTAGGGAAGATCTGGCCTTCCTACAAGGGAAGGCCAGGGAATTTTCTTCAGAGCAGACCAGA"
    "GCCAACAGCCCCACCAGAAGAGAGCTTCAGGTCTGGGGTAGAGACAACAACTCCCCCTCAGAAGCAGGAGCCGATAGACAAGGAAC"
    "TGTATCCTTTAACTTCCCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTAT"
    "TAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGT"
    "TTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACC"
    "TGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTTCCCATTAGCCCTATTGAGACTGTACCAGTAA"
    "AATTAAAGCCAGGAATGGATGGCCCAAAAGTTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACA"
    "GAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAG"
    "TACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCTGGGAAGTTCAATTAGGAATACCACATC"
    "CCGCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTGGGTGATGCATATTTTTCAGTTCCCTTAGATGAAGACTTCAGG"
    "AAGTATACTGCATTTACCATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAGGGATGGAA"
    "AGGATCACCAGCAATATTCCAAAGTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAAT"
    "ACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAGCTGAGACAACATCTGTTGAGG"
    "TGGGGACTTACCACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTTGGATGGGTTATGAACTCCATCCTGATAAATGGAC"
    "AGTACAGCCTATAGTGCTGCCAGAAAAAGACAGCTGGACTGTCAATGACATACAGAAGTTAGTGGGGAAATTGAATTGGGCAAGTC"
    "AGATTTACCCAGGGATTAAAGTAAGGCAATTATGTAAACTCCTTAGAGGAACCAAAGCACTAACAGAAGTAATACCACTAACAGAA"
    "GAAGCAGAGCTAGAACTGGCAGAAAACAGAGAGATTCTAAAAGAACCAGTACATGGAGTGTATTATGACCCATCAAAAGACTTAAT"
    "AGCAGAAATACAGAAGCAGGGGCAAGGCCAATGGACATATCAAATTTATCAAGAGCCATTTAAAAATCTGAAAACAGGAAAATATG"
    "CAAGAATGAGGGGTGCCCACACTAATGATGTAAAACAATTAACAGAGGCAGTGCAAAAAATAACCACAGAAAGCATAGTAATATGG"
    "GGAAAGACTCCTAAATTTAAACTGCCCATACAAAAGGAAACATGGGAAACATGGTGGACAGAGTATTGGCAAGCCACCTGGATTCC"
    "TGAGTGGGAGTTTGTTAATACCCCTCCCTTAGTGAAATTATGGTACCAGTTAGAGAAAGAACCCATAGTAGGAGCAGAAACCTTCT"
    "ATGTAGATGGGGCAGCTAACAGGGAGACTAAATTAGGAAAAGCAGGATATGTTACTAATAGAGGAAGACAAAAAGTTGTCACCCTA"
    "ACTGACACAACAAATCAGAAGACTGAGTTACAAGCAATTTATCTAGCTTTGCAGGATTCGGGATTAGAAGTAAACATAGTAACAGA"
    "CTCACAATATGCATTAGGAATCATTCAAGCACAACCAGATCAAAGTGAATCAGAGTTAGTCAATCAAATAATAGAGCAGTTAATAA"
    "AAAAGGAAAAGGTCTATCTGGCATGGGTACCAGCACACAAAGGAATTGGAGGAAATGAACAAGTAGATAAATTAGTCAGTGCTGGA"
    "ATCAGGAAAGTACTATTTTTAGATGGAATAGATAAGGCCCAAGATGAACATGAGAAATATCACAGTAATTGGAGAGCAATGGCTAG"
    "TGATTTTAACCTGCCACCTGTAGTAGCAAAAGAAATAGTAGCCAGCTGTGATAAATGTCAGCTAAAAGGAGAAGCCATGCATGGAC"
    "AAGTAGACTGTAGTCCAGGAATATGGCAACTAGATTGTACACATTTAGAAGGAAAAGTTATCCTGGTAGCAGTTCATGTAGCCAGT"
    "GGATATATAGAAGCAGAAGTTATTCCAGCAGAAACAGGGCAGGAAACAGCATATTTTCTTTTAAAATTAGCAGGAAGATGGCCAGT"
    "AAAAACAATACATACTGACAATGGCAGCAATTTCACCGGTGCTACGGTTAGGGCCGCCTGTTGGTGGGCGGGAATCAAGCAGGAAT"
    "TTGGAATTCCCTACAATCCCCAAAGTCAAGGAGTAGTAGAATCTATGAATAAAGAATTAAAGAAAATTATAGGACAGGTAAGAGAT"
    "CAGGCTGAACATCTTAAGACAGCAGTACAAATGGCAGTATTCATCCACAATTTTAAAAGAAAAGGGGGGATTGGGGGGTACAGTGC"
    "AGGGGAAAGAATAGTAGACATAATAGCAACAGACATACAAACTAAAGAATTACAAAAACAAATTACAAAAATTCAAAATTTTCGGG"
    "TTTATTACAGGGACAGCAGAAATCCACTTTGGAAAGGACCAGCAAAGCTCCTCTGGAAAGGTGAAGGGGCAGTAGTAATACAAGAT"
    "AATAGTGACATAAAAGTAGTGCCAAGAAGAAAAGCAAAGATCATTAGGGATTATGGAAAACAGATGGCAGGTGATGATTGTGTGGC"
    "AAGTAGACAGGATGAGGATTAGAACATGGAAAAGTTTAGTAAAACACCATATGTATGTTTCAGGGAAAGCTAGGGGATGGTTTTAT"
    "AGACATCACTATGAAAGCCCTCATCCAAGAATAAGTTCAGAAGTACACATCCCACTAGGGGATGCTAGATTGGTAATAACAACATA"
    "TTGGGGTCTGCATACAGGAGAAAGAGACTGGCATTTGGGTCAGGGAGTCTCCATAGAATGGAGGAAAAAGAGATATAGCACACAAG"
    "TAGACCCTGAACTAGCAGACCAACTAATTCATCTGTATTACTTTGACTGTTTTTCAGACTCTGCTATAAGAAAGGCCTTATTAGGA"
    "CACATAGTTAGCCCTAGGTGTGAATATCAAGCAGGACATAACAAGGTAGGATCTCTACAATACTTGGCACTAGCAGCATTAATAAC"
    "ACCAAAAAAGATAAAGCCACCTTTGCCTAGTGTTACGAAACTGACAGAGGATAGATGGAACAAGCCCCAGAAGACCAAGGGCCACA"
    "GAGGGAGCCACACAATGAATGGACACTAGAGCTTTTAGAGGAGCTTAAGAATGAAGCTGTTAGACATTTTCCTAGGATTTGGCTCC"
    "ATGGCTTAGGGCAACATATCTATGAAACTTATGGGGATACTTGGGCAGGAGTGGAAGCCATAATAAGAATTCTGCAACAACTGCTG"
    "TTTATCCATTTTCAGAATTGGGTGTCGACATAGCAGAATAGGCGTTACTCGACAGAGGAGAGCAAGAAATGGAGCCAGTAGATCCT"
    "AGACTAGAGCCCTGGAAGCATCCAGGAAGTCAGCCTAAAACTGCTTGTACCAATTGCTATTGTAAAAAGTGTTGCTTTCATTGCCA"
    "AGTTTGTTTCATAACAAAAGCCTTAGGCATCTCCTATGGCAGGAAGAAGCGGAGACAGCGACGAAGAGCTCATCAGAACAGTCAGA"
    "CTCATCAAGCTTCTCTATCAAAGCAGTAAGTAGTACATGTAACGCAACCTATACCAATAGTAGCAATAGTAGCATTAGTAGTAGCA"
    "ATAATAATAGCAATAGTTGTGTGGTCCATAGTAATCATAGAATATAGGAAAATATTAAGACAAAGAAAAATAGACAGGTTAATTGA"
    "TAGACTAATAGAAAGAGCAGAAGACAGTGGCAATGAGAGTGAAGGAGAAATATCAGCACTTGTGGAGATGGGGGTGGAGATGGGGC"
    "ACCATGCTCCTTGGGATGTTGATGATCTGTAGTGCTACAGAAAAATTGTGGGTCACAGTCTATTATGGGGTACCTGTGTGGAAGGA"
    "AGCAACCACCACTCTATTTTGTGCATCAGATGCTAAAGCATATGATACAGAGGTACATAATGTTTGGGCCACACATGCCTGTGTAC"
    "CCACAGACCCCAACCCACAAGAAGTAGTATTGGTAAATGTGACAGAAAATTTTAACATGTGGAAAAATGACATGGTAGAACAGATG"
    "CATGAGGATATAATCAGTTTATGGGATCAAAGCCTAAAGCCATGTGTAAAATTAACCCCACTCTGTGTTAGTTTAAAGTGCACTGA"
    "TTTGAAGAATGATACTAATACCAATAGTAGTAGCGGGAGAATGATAATGGAGAAAGGAGAGATAAAAAACTGCTCTTTCAATATCA"
    "GCACAAGCATAAGAGGTAAGGTGCAGAAAGAATATGCATTTTTTTATAAACTTGATATAATACCAATAGATAATGATACTACCAGC"
    "TATAAGTTGACAAGTTGTAACACCTCAGTCATTACACAGGCCTGTCCAAAGGTATCCTTTGAGCCAATTCCCATACATTATTGTGC"
    "CCCGGCTGGTTTTGCGATTCTAAAATGTAATAATAAGACGTTCAATGGAACAGGACCATGTACAAATGTCAGCACAGTACAATGTA"
    "CACATGGAATTAGGCCAGTAGTATCAACTCAACTGCTGTTAAATGGCAGTCTAGCAGAAGAAGAGGTAGTAATTAGATCTGTCAAT"
    "TTCACGGACAATGCTAAAACCATAATAGTACAGCTGAACACATCTGTAGAAATTAATTGTACAAGACCCAACAACAATACAAGAAA"
    "AAGAATCCGTATCCAGAGAGGACCAGGGAGAGCATTTGTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGTAACATTA"
    "GTAGAGCAAAATGGAATAACACTTTAAAACAGATAGCTAGCAAATTAAGAGAACAATTTGGAAATAATAAAACAATAATCTTTAAG"
    "CAATCCTCAGGAGGGGACCCAGAAATTGTAACGCACAGTTTTAATTGTGGAGGGGAATTTTTCTACTGTAATTCAACACAACTGTT"
    "TAATAGTACTTGGTTTAATAGTACTTGGAGTACTGAAGGGTCAAATAACACTGAAGGAAGTGACACAATCACCCTCCCATGCAGAA"
    "TAAAACAAATTATAAACATGTGGCAGAAAGTAGGAAAAGCAATGTATGCCCCTCCCATCAGTGGACAAATTAGATGTTCATCAAAT"
    "ATTACAGGGCTGCTATTAACAAGAGATGGTGGTAATAGCAACAATGAGTCCGAGATCTTCAGACCTGGAGGAGGAGATATGAGGGA"
    "CAATTGGAGAAGTGAATTATATAAATATAAAGTAGTAAAAATTGAACCATTAGGAGTAGCACCCACCAAGGCAAAGAGAAGAGTGG"
    "TGCAGAGAGAAAAAAGAGCAGTGGGAATAGGAGCTTTGTTCCTTGGGTTCTTGGGAGCAGCAGGAAGCACTATGGGCGCAGCCTCA"
    "ATGACGCTGACGGTACAGGCCAGACAATTATTGTCTGGTATAGTGCAGCAGCAGAACAATTTGCTGAGGGCTATTGAGGCGCAACA"
    "GCATCTGTTGCAACTCACAGTCTGGGGCATCAAGCAGCTCCAGGCAAGAATCCTGGCTGTGGAAAGATACCTAAAGGATCAACAGC"
    "TCCTGGGGATTTGGGGTTGCTCTGGAAAACTCATTTGCACCACTGCTGTGCCTTGGAATGCTAGTTGGAGTAATAAATCTCTGGAA"
    "CAGATTTGGAATCACACGACCTGGATGGAGTGGGACAGAGAAATTAACAATTACACAAGCTTAATACACTCCTTAATTGAAGAATC"
    "GCAAAACCAGCAAGAAAAGAATGAACAAGAATTATTGGAATTAGATAAATGGGCAAGTTTGTGGAATTGGTTTAACATAACAAATT"
    "GGCTGTGGTATATAAAATTATTCATAATGATAGTAGGAGGCTTGGTAGGTTTAAGAATAGTTTTTGCTGTACTTTCTATAGTGAAT"
    "AGAGTTAGGCAGGGATATTCACCATTATCGTTTCAGACCCACCTCCCAACCCCGAGGGGACCCGACAGGCCCGAAGGAATAGAAGA"
    "AGAAGGTGGAGAGAGAGACAGAGACAGATCCATTCGATTAGTGAACGGATCCTTGGCACTTATCTGGGACGATCTGCGGAGCCTGT"
    "GCCTCTTCAGCTACCACCGCTTGAGAGACTTACTCTTGATTGTAACGAGGATTGTGGAACTTCTGGGACGCAGGGGGTGGGAAGCC"
    "CTCAAATATTGGTGGAATCTCCTACAGTATTGGAGTCAGGAACTAAAGAATAGTGCTGTTAGCTTGCTCAATGCCACAGCCATAGC"
    "AGTAGCTGAGGGGACAGATAGGGTTATAGAAGTAGTACAAGGAGCTTGTAGAGCTATTCGCCACATACCTAGAAGAATAAGACAGG"
    "GCTTGGAAAGGATTTTGCTATAAGATGGGTGGCAAGTGGTCAAAAAGTAGTGTGATTGGATGGCCTACTGTAAGGGAAAGAATGAG"
    "ACGAGCTGAGCCAGCAGCAGATAGGGTGGGAGCAGCATCTCGAGACCTGGAAAAACATGGAGCAATCACAAGTAGCAATACAGCAG"
    "CTACCAATGCTGCTTGTGCCTGGCTAGAAGCACAAGAGGAGGAGGAGGTGGGTTTTCCAGTCACACCTCAGGTACCTTTAAGACCA"
    "ATGACTTACAAGGCAGCTGTAGATCTTAGCCACTTTTTAAAAGAAAAGGGGGGACTGGAAGGGCTAATTCACTCCCAAAGAAGACA"
    "AGATATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACCAGGGCCAGGGGTCAGATATC"
    "CACTGACCTTTGGATGGTGCTACAAGCTAGTACCAGTTGAGCCAGATAAGATAGAAGAGGCCAATAAAGGAGAGAACACCAGCTTG"
    "TTACACCCTGTGAGCCTGCATGGGATGGATGACCCGGAGAGAGAAGTGTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCA"
    "CGTGGCCCGAGAGCTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAGGGACTTTCCGCTGGGGACTTTCCAG"
    "GGAGGCGTGGCCTGGGCGGGACTGGGGAGTGGCGAGCCCTCAGATCCTGCATATAAGCAGCTGCTTTTTGCCTGTACTGGGTCTCT"
    "CTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAGCCTCAATAAAGCTTGCCTTGAGTGC"
    "TTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGC"
    "A";

const std::vector<int> AminoAcidCaller::nnrti = {100, 101, 103, 106, 138, 179,
                                                 181, 190, 190, 227, 230};
const std::vector<int> AminoAcidCaller::nnrtiSurveillance = {41,  65,  67,  69,  70,  74,  75, 77,
                                                             115, 116, 151, 184, 210, 215, 219};
const std::vector<int> AminoAcidCaller::nrtiSurveillance = {100, 101, 103, 106, 179,
                                                            181, 188, 190, 225, 230};
const std::vector<int> AminoAcidCaller::nrti = {184, 65,  70,  74,  115, 41, 67,
                                                70,  210, 215, 219, 69,  151};
const std::vector<int> AminoAcidCaller::pi = {24, 32, 46, 47, 50, 54, 76, 82, 84, 88, 90};
const std::vector<int> AminoAcidCaller::piSurveillance = {23, 24, 30, 32, 46, 47, 48, 50, 53,
                                                          54, 73, 76, 82, 83, 84, 85, 88, 90};
const std::vector<int> AminoAcidCaller::ini = {66, 92, 138, 140, 143, 147, 148, 155};
}
}  // ::PacBio::Juliet
