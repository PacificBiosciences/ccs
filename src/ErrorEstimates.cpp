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

#include <pacbio/juliet/ErrorEstimates.h>

namespace PacBio {
namespace Juliet {

ErrorEstimates::ErrorEstimates(const std::string& s) { SetFromModel(ErrorModelFromString(s)); }

ErrorEstimates::ErrorEstimates(const ErrorModel& e) { SetFromModel(e); }

void ErrorEstimates::SetFromModel(const ErrorModel& e)
{
    switch (e) {
        case ErrorModel::SP1C1_RQ99:
            match = 0.9930786;
            substitution = 0.0007421148 / 3.0;  // 0.0006101725 + 3*4.398076e-05
            deletion = 0.006179274;             // 0.003515625 + 3*0.0008878829
            insertion = 0;
            break;
        case ErrorModel::SP1C1_RQ95:
            match = 0.9877258;
            substitution = 0.00216356 / 3.0;  // 0.001664215 + 3*0.0001664483
            deletion = 0.01011063;            // 0.00646245 + 3*0.001216059
            insertion = 0;
            break;
        default:
            throw std::runtime_error("Unknown error model");
    }
}
}
}