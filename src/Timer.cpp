// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
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

// Author: Lance Hepler, Armin Töpfer

#include <chrono>
#include <sstream>

#include <pacbio/util/Timer.h>

using std::chrono::steady_clock;
using std::chrono::milliseconds;

namespace PacBio {
namespace Util {

Timer::Timer() { Restart(); }

float Timer::ElapsedMilliseconds() const
{
    auto tock = steady_clock::now();
    return std::chrono::duration_cast<milliseconds>(tock - tick).count();
}

float Timer::ElapsedSeconds() const { return ElapsedMilliseconds() / 1000; }

std::string Timer::ElapsedTime() const
{
    auto tock = steady_clock::now();
    auto t = std::chrono::duration_cast<std::chrono::nanoseconds>(tock - tick).count();

    auto d = t / 10000 / 1000 / 1000 / 60 / 60 / 24;
    auto h = (t / 1000 / 1000 / 1000 / 60 / 60) % 24;
    auto m = (t / 1000 / 1000 / 1000 / 60) % 60;
    auto s = (t / 1000 / 1000 / 1000) % 60;
    auto ms = (t / 1000 / 1000) % 1000;
    auto us = (t / 1000) % 1000;
    auto ns = t % 1000;
    std::stringstream ss;
    if (d > 0) ss << d << "d ";
    if (h > 0) ss << h << "h ";
    if (m > 0) ss << m << "m ";
    if (s > 0) ss << s << "s ";
    if (ms > 0) ss << ms << "ms ";
    if (us > 0) ss << us << "us ";
    if (ns > 0) ss << ns << "ns ";
    return ss.str();
}

void Timer::Restart() { tick = steady_clock::now(); }

}  // namespace Util
}  // namespace PacBio
