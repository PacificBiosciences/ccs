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

#include <array>

#include "RandomDNA.h"

using namespace std;

constexpr double pulseWidthDist[50] = {
    0.0,         0.0,         67.58137063, 72.67947116, 71.90971706, 68.52503331, 63.93130097,
    58.79290293, 53.5441654,  48.47127299, 43.69099241, 39.30492681, 35.29555429, 31.69712306,
    28.45523289, 25.56406112, 22.96958546, 20.66754424, 18.61712377, 16.78288686, 15.15159641,
    13.69738276, 12.39031514, 11.22761988, 10.1905659,  9.25254787,  8.41958363,  7.66550808,
    6.99112275,  6.37780398,  5.82998074,  5.33196501,  4.88806118,  4.47817299,  4.1106286,
    3.78052439,  3.47567456,  3.20017765,  2.94609885,  2.72085145,  2.51190674,  2.32144746,
    2.14614522,  1.9881507,   1.84008447,  1.70632031,  1.58504037,  1.47055015,  1.36683485,
    1.27025397};

string RandomDNA(const size_t n, mt19937* const gen)
{
    static constexpr const std::array<char, 4> bases{{'A', 'C', 'G', 'T'}};

    string result(n, 'A');
    uniform_int_distribution<size_t> rand(0, 3);

    for (size_t i = 0; i < n; ++i)
        result[i] = bases[rand(*gen)];

    return result;
}

std::vector<uint8_t> RandomPW(const size_t n, mt19937* const gen)
{
    std::vector<uint8_t> result(n);
    discrete_distribution<uint8_t> rand(begin(pulseWidthDist), end(pulseWidthDist));
    for (size_t i = 0; i < n; ++i)
        result[i] = rand(*gen);

    return result;
}
