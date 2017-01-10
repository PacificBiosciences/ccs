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

#include <algorithm>
#include <exception>
#include <locale>
#include <string>

namespace PacBio {
namespace Juliet {

enum class ErrorModel : uint8_t
{
    SP1C1_RQ95 = 0,
    SP1C1_RQ99
};

inline ErrorModel ErrorModelFromString(const std::string& input)
{
    ErrorModel e;
    std::string s = input;
    std::transform(s.begin(), s.end(), s.begin(), ::toupper);
    if (s == "SP1C1_RQ95")
        e = ErrorModel::SP1C1_RQ95;
    else if (s == "SP1C1_RQ99")
        e = ErrorModel::SP1C1_RQ99;
    else
        throw std::runtime_error("Unknown error model string");
    return e;
}

/// Contains CCS error estimates
class ErrorEstimates
{
public:
    ErrorEstimates(const std::string& s);
    ErrorEstimates(const ErrorModel& e);

    double match;
    double substitution;
    double deletion;
    double insertion;

private:
    void SetFromModel(const ErrorModel& e);
};
}
}  // ::PacBio::Juliet