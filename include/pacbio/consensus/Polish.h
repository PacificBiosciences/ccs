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

#include <tuple>
#include <vector>

#include <pacbio/consensus/Mutation.h>
#include <pacbio/consensus/PolishResult.h>

namespace PacBio {
namespace Consensus {

// forward declaration
class Integrator;

struct PolishConfig
{
    size_t MaximumIterations;
    size_t MutationSeparation;
    size_t MutationNeighborhood;

    PolishConfig(size_t iterations = 40, size_t separation = 10, size_t neighborhood = 20);
};

struct RepeatConfig
{
    size_t MaximumRepeatSize;
    size_t MinimumElementCount;
    size_t MaximumIterations;

    RepeatConfig(size_t repeatSize = 3, size_t elementCount = 3, size_t iterations = 40);
};

/// Given an AbstractIntegrator and a PolishConfig,
/// Given an Integrator and a PolishConfig,
/// iteratively polish the template,
/// and return meta information about the procedure.
///
/// The template will be polished within the Integrator.
PolishResult Polish(Integrator* ai, const PolishConfig& cfg);

PolishResult PolishRepeats(AbstractIntegrator* ai, const RepeatConfig& cfg);

/// Struct that contains vectors for the base-wise individual and compound QVs.
struct QualityValues
{
    std::vector<int> Qualities;
    std::vector<int> DeletionQVs;
    std::vector<int> InsertionQVs;
    std::vector<int> SubstitutionQVs;
};

/// Generates phred qualities of the current template.
std::vector<int> ConsensusQualities(Integrator& ai);

/// Generates individual and compound phred qualities of the current template.
QualityValues ConsensusQVs(Integrator& ai);

/// Returns a list of all possible mutations that can be applied to the template
/// of the provided integrator.
std::vector<Mutation> Mutations(const Integrator& ai);

/// Returns a list of all possible repeat mutations of the template
/// of the provided integrator
std::vector<Mutation> RepeatMutations(const AbstractIntegrator& ai, const RepeatConfig& cfg);

}  // namespace Consensus
}  // namespace PacBio
