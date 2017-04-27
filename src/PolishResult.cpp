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

#include <pacbio/consensus/PolishResult.h>

namespace PacBio {
namespace Consensus {

PolishResult operator+(const PolishResult& lhs, const PolishResult& rhs)
{
    PolishResult result;
    result.hasConverged = lhs.hasConverged && rhs.hasConverged;
    result.mutationsTested = lhs.mutationsTested + rhs.mutationsTested;
    result.mutationsApplied = lhs.mutationsApplied + rhs.mutationsApplied;
    result.maxAlphaPopulated.insert(result.maxAlphaPopulated.end(), lhs.maxAlphaPopulated.begin(),
                                    lhs.maxAlphaPopulated.end());
    result.maxBetaPopulated.insert(result.maxBetaPopulated.end(), lhs.maxBetaPopulated.begin(),
                                   lhs.maxBetaPopulated.end());
    result.maxNumFlipFlops.insert(result.maxNumFlipFlops.end(), lhs.maxNumFlipFlops.begin(),
                                  lhs.maxNumFlipFlops.end());
    result.maxAlphaPopulated.insert(result.maxAlphaPopulated.end(), rhs.maxAlphaPopulated.begin(),
                                    rhs.maxAlphaPopulated.end());
    result.maxBetaPopulated.insert(result.maxBetaPopulated.end(), rhs.maxBetaPopulated.begin(),
                                   rhs.maxBetaPopulated.end());
    result.maxNumFlipFlops.insert(result.maxNumFlipFlops.end(), rhs.maxNumFlipFlops.begin(),
                                  rhs.maxNumFlipFlops.end());
    return result;
}
}
}  // ::PacBio::Consensus
