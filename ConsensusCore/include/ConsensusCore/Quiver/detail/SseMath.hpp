// Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.
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

// Author: David Alexander

#pragma once

#include <xmmintrin.h>
#include <limits>

#include <ConsensusCore/Quiver/detail/sse_mathfun.h>


// todo: turn these into inline functions
#define ADD4(a, b) _mm_add_ps((a), (b))

#define AFFINE4(offset, slope, dataptr)                 \
  (_mm_add_ps(_mm_set_ps1(offset),                      \
              _mm_mul_ps(_mm_set_ps1(slope),            \
                         _mm_loadu_ps((dataptr)))))

#define MUX4(mask, a, b) (_mm_or_ps(_mm_and_ps((mask), (a)), _mm_andnot_ps((mask), (b))))

#define MAX4(a, b) _mm_max_ps((a), (b))


namespace ConsensusCore {
namespace detail {
    //
    // Log-space arithmetic
    //
    static const __m128 ones    = _mm_set_ps1(1.0f);

    inline __m128 logAdd4(__m128 aa, __m128 bb)
    {
        __m128 max = _mm_max_ps(aa, bb);
        __m128 min = _mm_min_ps(aa, bb);
        __m128 diff = _mm_sub_ps(min, max);
        return _mm_add_ps(max, log_ps(_mm_add_ps(ones, exp_ps(diff))));

        // Precision of this is currently too low -- causes problems w/ Edna
        // need a more accurate impl.
        // return logAddApprox_ps(aa, bb);
    }

    inline float logAdd(float a, float b)
    {
        __m128 aa = _mm_set_ps1(a);
        __m128 bb = _mm_set_ps1(b);
        __m128 acc = logAdd4(aa, bb);

        // Use an aligned buf to extract
        ALIGN16_BEG float buf[4] ALIGN16_END;
        _mm_store_ps(&buf[0], acc);
        return buf[0];
    }
}}
