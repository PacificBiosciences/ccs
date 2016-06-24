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

#include "Mutations.h"

std::vector<PacBio::Consensus::Mutation> Mutations(const std::string& tpl, const size_t start,
                                                   const size_t end)
{
    using namespace PacBio::Consensus;

    constexpr auto bases = "ACGT";

    std::vector<Mutation> result;

    for (size_t i = start; i < end; ++i) {
        for (size_t j = 0; j < 4; ++j)
            result.push_back(Mutation(MutationType::INSERTION, i, bases[j]));

        result.push_back(Mutation(MutationType::DELETION, i));

        for (size_t j = 0; j < 4; ++j)
            if (bases[j] != tpl[i])
                result.push_back(Mutation(MutationType::SUBSTITUTION, i, bases[j]));
    }

    for (size_t j = 0; j < 4; ++j)
        result.push_back(Mutation(MutationType::INSERTION, tpl.length(), bases[j]));

    return result;
}

std::vector<PacBio::Consensus::Mutation> Mutations(const std::string& tpl)
{
    return Mutations(tpl, 0, tpl.length());
}
