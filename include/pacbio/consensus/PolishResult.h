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

#pragma once

#include <vector>

#include <stddef.h>

#include <boost/optional.hpp>

#include <pacbio/consensus/Mutation.h>

namespace PacBio {
namespace Consensus {

struct DiploidSite
{
    MutationType mutType;
    std::vector<char> mutants;
    int64_t pos;
    boost::optional<double> pvalue;

    DiploidSite(const MutationType mutType_, const std::vector<char>& mutants_, const int64_t pos_,
                const boost::optional<double> pvalue_ = boost::none)
        : mutType{mutType_}, mutants{mutants_}, pos{pos_}, pvalue{pvalue_}
    {
    }

    bool operator==(const DiploidSite& rhs) const
    {
        return std::tie(mutType, mutants, pos, pvalue) ==
               std::tie(rhs.mutType, rhs.mutants, rhs.pos, rhs.pvalue);
    }
};

/// This struct contains the results of Integrator::Polish()
struct PolishResult
{
    // Did Polish() converge?
    bool hasConverged = false;
    // How many mutations have been tested?
    size_t mutationsTested = 0;
    // How many mutations have been actually applied?
    size_t mutationsApplied = 0;

    // For each iteration in Polish(), get the max of all Evaluators to
    // diagnose the worst performing one.
    //
    // Maximal ratio of populated alpha cells
    std::vector<float> maxAlphaPopulated;
    // Maximal ratio of populated beta cells
    std::vector<float> maxBetaPopulated;
    // Maximal number of flip flop events
    std::vector<int> maxNumFlipFlops;

    // Diploid results
    // The vector is sorted according to the standard
    // unanimity Mutation class criterion
    std::vector<DiploidSite> diploidSites;
};

PolishResult operator+(const PolishResult& lhs, const PolishResult& rhs);
}
}  // ::PacBio::Consensus