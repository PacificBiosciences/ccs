// Copyright (c) 2011-2017, Pacific Biosciences of California, Inc.
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

#include <cmath>
#include <functional>

#include <pacbio/consensus/ModelConfig.h>

namespace PacBio {
namespace Consensus {
namespace {

inline double CounterWeight(std::function<double(size_t, MoveType)>&& lgPrTransition,
                            std::function<double(size_t, MoveType)>&& lgPrEmission,
                            const size_t nContexts)
{
    double meanPrEm = 0.0;
    for (size_t ctx = 0; ctx < nContexts; ++ctx) {
        const double pr_M = lgPrTransition(ctx, MoveType::MATCH);
        const double pr_B = lgPrTransition(ctx, MoveType::BRANCH);
        const double pr_S = lgPrTransition(ctx, MoveType::STICK);
        const double pr_D = lgPrTransition(ctx, MoveType::DELETION);

        const double lgPr_EM = lgPrEmission(ctx, MoveType::MATCH);
        const double lgPr_EB = lgPrEmission(ctx, MoveType::BRANCH);
        const double lgPr_ES = lgPrEmission(ctx, MoveType::STICK);
        const double lgPr_ED = 0.0;  // nothing to emit

        const double E_MD = lgPr_EM * pr_M / (pr_M + pr_D) + lgPr_ED * pr_D / (pr_M + pr_D);
        const double E_BS = lgPr_EB * pr_B / (pr_B + pr_S) + lgPr_ES * pr_S / (pr_B + pr_S);
        const double E_I = E_BS * (pr_B + pr_S) / (pr_M + pr_D);

        meanPrEm += std::exp(E_I + E_MD);
    }
    meanPrEm /= nContexts;
    return 1.0 / meanPrEm;
}
}
}
}
