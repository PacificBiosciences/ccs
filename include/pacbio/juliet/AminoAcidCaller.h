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

#include <array>
#include <cmath>
#include <exception>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <list>
#include <memory>
#include <numeric>
#include <sstream>
#include <unordered_map>
#include <vector>

#include <pacbio/data/MSA.h>
#include <pbcopper/json/JSON.h>

namespace PacBio {
namespace Juliet {
namespace {
struct VariantGene
{
    std::string geneName;

    struct VariantPosition
    {
        std::string refCodon;
        char refAminoAcid;
        std::vector<JSON::Json> msa;
        int coverage;

        struct VariantCodon
        {
            std::string codon;
            double frequency;
            double pValue;
            std::string knownDRM;
        };
        std::map<char, std::vector<VariantCodon>> aminoAcidToCodons;
    };

    std::map<int, VariantPosition> relPositionToVariant;

    JSON::Json ToJson() const
    {
        using namespace JSON;
        Json root;
        root["name"] = geneName;
        std::vector<Json> positions;
        for (const auto& pos_variant : relPositionToVariant) {
            Json jVarPos;
            jVarPos["ref_position"] = pos_variant.first;
            jVarPos["ref_codon"] = pos_variant.second.refCodon;
            jVarPos["coverage"] = pos_variant.second.coverage;
            jVarPos["ref_amino_acid"] = std::string(1, pos_variant.second.refAminoAcid);

            if (pos_variant.second.aminoAcidToCodons.empty()) continue;
            std::vector<Json> jVarAAs;
            for (const auto& aa_varCodon : pos_variant.second.aminoAcidToCodons) {
                Json jVarAA;
                jVarAA["amino_acid"] = std::string(1, aa_varCodon.first);
                std::vector<Json> jCodons;

                if (aa_varCodon.second.empty()) continue;
                for (const auto& codon : aa_varCodon.second) {
                    Json jCodon;
                    jCodon["codon"] = codon.codon;
                    jCodon["frequency"] = codon.frequency;
                    jCodon["pValue"] = codon.pValue;
                    jCodon["known_drm"] = codon.knownDRM;
                    jCodons.push_back(jCodon);
                }
                jVarAA["variant_codons"] = jCodons;
                jVarAAs.push_back(jVarAA);
            }
            jVarPos["variant_amino_acids"] = jVarAAs;
            jVarPos["msa"] = pos_variant.second.msa;
            positions.push_back(jVarPos);
        }
        if (!positions.empty()) root["variant_positions"] = positions;
        return root;
    }
};
}

/// Given a MSA and p-values for each nucleotide of each position,
/// generate machine-interpretable and human-readable output about mutated
/// amino acids.
class AminoAcidCaller
{
public:
    AminoAcidCaller(const std::vector<Data::ArrayRead>& reads);

public:
    /// Generate JSON output of variant amino acids
    JSON::Json JSON();

public:
    /// Generate HTML output of variant amino acids
    static void HTML(std::ostream& out, const JSON::Json& j, bool onlyKnownDRMs, bool details);

public:
    std::unique_ptr<Data::MSA> msa_;

private:
    void GenerateMSA(const std::vector<Data::ArrayRead>& reads);
    void CallVariants(const std::vector<Data::ArrayRead>& reads);

private:
    int beginPos_ = std::numeric_limits<int>::max();
    int endPos_ = 0;
    std::vector<std::vector<char>> matrix_;
    std::vector<VariantGene> variantGenes_;

    static const std::unordered_map<std::string, char> codonToAmino_;
    static const std::string ref_;
    static const std::unordered_map<int, std::string> resistentCodon_;

    static const std::vector<int> nnrti;
    static const std::vector<int> nnrtiSurveillance;
    static const std::vector<int> nrti;
    static const std::vector<int> nrtiSurveillance;
    static const std::vector<int> pi;
    static const std::vector<int> piSurveillance;
    static const std::vector<int> ini;
};
}
}  // ::PacBio::Juliet