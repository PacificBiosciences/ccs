/*
 * This is a port of some functions from the Yepp library. Original copywrite below.
 *
 *                          Yeppp! library header
 *
 * This file is part of Yeppp! library and licensed under the New BSD license.
 *
 * Copyright (C) 2010-2012 Marat Dukhan
 * Copyright (C) 2012-2013 Georgia Institute of Technology
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Georgia Institute of Technology nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#pragma once

#include <stdbool.h>
#include <stdint.h>
#include <x86intrin.h>

#define YEP_INLINE inline
#define YEP_NATIVE_FUNCTION
#define YEP_LIKELY(x) (__builtin_expect(!!(x), 1))
#define YEP_UNLIKELY(x) (__builtin_expect(!!(x), 0))
#define YEP_PROCESSOR_SUPPORTS_FMA_EXTENSION 1

// ND - I think this is true for intel after 2013, not going to define it for now
// though
//#define YEP_PROCESSOR_SUPPORTS_FMA_EXTENSION

typedef double                         Yep64f;
typedef uint64_t                       Yep64u;
typedef uint32_t                       Yep32u;
typedef int32_t                        Yep32s;
typedef bool               YepBoolean;
typedef struct {
    Yep64f high;
    Yep64f low;
} Yep64df;


YEP_NATIVE_FUNCTION static YEP_INLINE Yep32s yepBuiltin_Clamp_32s32s32s_32s(Yep32s x, Yep32s xMin, Yep32s xMax) {
    return (x < xMin) ? xMin : (x > xMax) ? xMax : x;
}

YEP_NATIVE_FUNCTION static YEP_INLINE YepBoolean yepBuiltin_IsNaN_64f(Yep64f n) {
    return !(n == n);
}

YEP_NATIVE_FUNCTION static YEP_INLINE Yep64u yepBuiltin_GetLowPart_64u_32u(Yep64u n) {
    return Yep32u(n);
}

YEP_NATIVE_FUNCTION static YEP_INLINE Yep64u yepBuiltin_CombineParts_32u32u_64u(Yep32u hi, Yep32u lo) {
    return (Yep64u(hi) << 32) | Yep64u(lo);
}

YEP_NATIVE_FUNCTION static YEP_INLINE Yep64u yepBuiltin_Cast_64f_64u(Yep64f x) {
#if defined(YEP_NVIDIA_COMPILER)
    return __double_as_longlong(x);
#elif defined(YEP_INTEL_COMPILER)
    return _castf64_u64(x);
#else
    union {
        Yep64f float64;
        Yep64u word64;
    } float64_word64;
    float64_word64.float64 = x;
    return float64_word64.word64;
#endif
}

YEP_NATIVE_FUNCTION static YEP_INLINE Yep64f yepBuiltin_Cast_64u_64f(Yep64u x) {
#if defined(YEP_NVIDIA_COMPILER)
    return __longlong_as_double(x);
#elif defined(YEP_INTEL_COMPILER)
    return _castu64_f64(x);
#else
    union {
        Yep64f float64;
        Yep64u word64;
    } float64_word64;
    float64_word64.word64 = x;
    return float64_word64.float64;
#endif
}

YEP_NATIVE_FUNCTION static YEP_INLINE Yep64f yepBuiltin_PositiveInfinity_64f() {
#if defined(YEP_GNU_COMPILER) || defined(YEP_CLANG_COMPILER) || defined(YEP_NVIDIA_COMPILER)
    return __builtin_inf();
#else
    const static Yep64f one = 1.0;
    const static Yep64f zero = 0.0;
    return one / zero;
#endif
}

#if defined(YEP_PROCESSOR_SUPPORTS_DOUBLE_PRECISION_FMA_INSTRUCTIONS)
YEP_NATIVE_FUNCTION static YEP_INLINE Yep64f yepBuiltin_Divide_64f64f64f_64f(Yep64f y, Yep64f c, Yep64f rcpC) {
    const Yep64f q = y * rcpC;
    const Yep64f r = yepBuiltin_FNMA_64f64f64f_64f(c, q, y);
    return yepBuiltin_FMA_64f64f64f_64f(r, rcpC, q);
}

YEP_NATIVE_FUNCTION static YEP_INLINE Yep64f yepBuiltin_MultiplyAdd_64f64f64f_64f(Yep64f a, Yep64f b, Yep64f c) {
    return yepBuiltin_FMA_64f64f64f_64f(a, b, c);
}

YEP_NATIVE_FUNCTION static YEP_INLINE Yep64f yepBuiltin_MultiplySubtract_64f64f64f_64f(Yep64f a, Yep64f b, Yep64f c) {
    return yepBuiltin_FMS_64f64f64f_64f(a, b, c);
}
#else
YEP_NATIVE_FUNCTION static YEP_INLINE Yep64f yepBuiltin_Divide_64f64f64f_64f(Yep64f y, Yep64f c, Yep64f rcpC) {
    return y / c;
}

YEP_NATIVE_FUNCTION static YEP_INLINE Yep64f yepBuiltin_MultiplyAdd_64f64f64f_64f(Yep64f a, Yep64f b, Yep64f c) {
    return __builtin_fma(a, b, c);
    //return a * b + c;
}

YEP_NATIVE_FUNCTION static YEP_INLINE Yep64f yepBuiltin_MultiplySubtract_64f64f64f_64f(Yep64f a, Yep64f b, Yep64f c) {
    return a * b - c;
}
#endif


YEP_NATIVE_FUNCTION static YEP_INLINE Yep64f yepBuiltin_Exp_64f_64f(Yep64f x) {
    
#if defined(YEP_COMPILER_SUPPORTS_HEXADECIMAL_FLOATING_POINT_CONSTANTS)
    const Yep64f magicBias = 0x1.8000000000000p+52;
    const Yep64f log2e     = 0x1.71547652B82FEp+0;
#if defined(YEP_PROCESSOR_SUPPORTS_FMA_EXTENSION)
    const Yep64df ln2 = { 0x1.62E42FEFA39EFp-1, 0x1.ABC9E3B39803Fp-56 };
#else
    const Yep64df ln2 = { 0x1.62E42FEFA3800p-1, 0x1.EF35793C76730p-45 };
#endif
    const Yep64f c2  = 0x1.0000000000005p-1;
    const Yep64f c3  = 0x1.5555555555540p-3;
    const Yep64f c4  = 0x1.5555555552115p-5;
    const Yep64f c5  = 0x1.11111111173CAp-7;
    const Yep64f c6  = 0x1.6C16C17F2BF99p-10;
    const Yep64f c7  = 0x1.A01A017EEB164p-13;
    const Yep64f c8  = 0x1.A019A6AC02A7Dp-16;
    const Yep64f c9  = 0x1.71DE71651CE7Ap-19;
    const Yep64f c10 = 0x1.28A284098D813p-22;
    const Yep64f c11 = 0x1.AE9043CA87A40p-26;
    
    const Yep64f zeroCutoff = -0x1.74910D52D3051p+9;
    const Yep64f infCutoff = 0x1.62E42FEFA39EFp+9;
#else
    const Yep64f magicBias = 6755399441055744.0;
    const Yep64f log2e     = 1.4426950408889634;
#if defined(YEP_PROCESSOR_SUPPORTS_FMA_EXTENSION)
    const Yep64df ln2 = { 0.6931471805599453, 2.3190468138462996e-17 };
#else
    const Yep64df ln2 = { 0.6931471805598903, 5.497923018708371e-14 };
#endif
    const Yep64f c2  = 0.5000000000000006;
    const Yep64f c3  = 0.16666666666666607;
    const Yep64f c4  = 0.04166666666657385;
    const Yep64f c5  = 0.008333333333377175;
    const Yep64f c6  = 0.0013888888932278352;
    const Yep64f c7  = 0.0001984126974695729;
    const Yep64f c8  = 2.4801504579877947e-5;
    const Yep64f c9  = 2.755738182142102e-6;
    const Yep64f c10 = 2.762627110160372e-7;
    const Yep64f c11 = 2.5062096212675488e-8;
    
    const Yep64f zeroCutoff = -745.1332191019411;
    const Yep64f infCutoff = 709.7827128933840;
#endif
    
    if YEP_UNLIKELY(yepBuiltin_IsNaN_64f(x)) {
        return x;
    } else {
        Yep64f t = x * log2e + magicBias;
        Yep32u e1 = yepBuiltin_GetLowPart_64u_32u(yepBuiltin_Cast_64f_64u(t)) << 20;
        Yep32u e2 = e1;
        e1 = yepBuiltin_Clamp_32s32s32s_32s(e1, -1022 << 20, 1023 << 20);
        e2 -= e1;
        const Yep64f s1 = yepBuiltin_Cast_64u_64f(yepBuiltin_CombineParts_32u32u_64u(e1 + 0x3FF00000u, 0u));
        const Yep64f s2 = yepBuiltin_Cast_64u_64f(yepBuiltin_CombineParts_32u32u_64u(e2 + 0x3FF00000u, 0u));
        t -= magicBias;
        const Yep64f rx = (x - t * ln2.high) - t * ln2.low;
        const Yep64f px = yepBuiltin_MultiplyAdd_64f64f64f_64f(rx,
                                                               yepBuiltin_MultiplyAdd_64f64f64f_64f(rx,
                                                                                                    yepBuiltin_MultiplyAdd_64f64f64f_64f(rx,
                                                                                                                                         yepBuiltin_MultiplyAdd_64f64f64f_64f(rx,
                                                                                                                                                                              yepBuiltin_MultiplyAdd_64f64f64f_64f(rx,
                                                                                                                                                                                                                   yepBuiltin_MultiplyAdd_64f64f64f_64f(rx,
                                                                                                                                                                                                                                                        yepBuiltin_MultiplyAdd_64f64f64f_64f(rx,
                                                                                                                                                                                                                                                                                             yepBuiltin_MultiplyAdd_64f64f64f_64f(rx,
                                                                                                                                                                                                                                                                                                                                  yepBuiltin_MultiplyAdd_64f64f64f_64f(rx, c11, c10),
                                                                                                                                                                                                                                                                                                                                  c9),
                                                                                                                                                                                                                                                                                             c8),
                                                                                                                                                                                                                                                        c7),
                                                                                                                                                                                                                   c6),
                                                                                                                                                                              c5),
                                                                                                                                         c4),
                                                                                                    c3),
                                                               c2);
        const Yep64f rf = yepBuiltin_MultiplyAdd_64f64f64f_64f(rx, rx * px, rx);
        Yep64f f = s2 * yepBuiltin_MultiplyAdd_64f64f64f_64f(s1, rf, s1);
        if YEP_UNLIKELY(x > infCutoff) {
            f = yepBuiltin_PositiveInfinity_64f();
        }
        if YEP_UNLIKELY(x < zeroCutoff) {
            f = 0.0;
        }
        return f;
    }
}

#endif
