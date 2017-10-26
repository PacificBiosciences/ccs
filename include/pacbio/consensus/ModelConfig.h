// Copyright (c) 2011-2016, Pacific Biosciences of California, Inc.
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

#pragma once

#include <array>
#include <cstdint>
#include <map>
#include <memory>
#include <ostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include <pacbio/data/Read.h>
#include <pacbio/data/internal/BaseEncoding.h>

namespace PacBio {

namespace Data {
struct MappedRead;
struct SNR;
}

namespace Consensus {

// fwd decl
class AbstractRecursor;
class AbstractTemplate;

using Data::detail::NCBI2na;
using Data::detail::NCBI4na;

// The allele representation in Unanimity
// currently employs the NCBI4na model, in
// order to account for diploid sites.
typedef NCBI4na AlleleRep;

struct TemplatePosition
{
    char Base;
    AlleleRep Idx;
    double Match;
    double Branch;
    double Stick;
    double Deletion;

public:
    // Constructor for backwards-compatibility
    constexpr TemplatePosition(const char base, const double match, const double branch,
                               const double stick, const double deletion)
        : Base{base}
        , Idx{AlleleRep::FromASCII(Base)}
        , Match{match}
        , Branch{branch}
        , Stick{stick}
        , Deletion{deletion}
    {
    }
};

std::ostream& operator<<(std::ostream&, const TemplatePosition&);

enum struct MoveType : uint8_t
{
    MATCH = 0,
    BRANCH = 1,
    STICK = 2,
    DELETION = 3  // never used for covariate
};

enum struct MomentType : uint8_t
{
    FIRST = 0,
    SECOND = 1
};

class ModelConfig
{
public:
    virtual ~ModelConfig() {}
    virtual std::unique_ptr<AbstractRecursor> CreateRecursor(const PacBio::Data::MappedRead& mr,
                                                             double scoreDiff) const = 0;
    virtual std::vector<TemplatePosition> Populate(const std::string& tpl) const = 0;
    virtual std::pair<Data::Read, std::vector<MoveType>> SimulateRead(
        std::default_random_engine* const rng, const std::string& tpl,
        const std::string& readname) const = 0;
    virtual double ExpectedLLForEmission(MoveType move, const AlleleRep& prev,
                                         const AlleleRep& curr, MomentType moment) const = 0;
};

}  // namespace Consensus
}  // namespace PacBio
