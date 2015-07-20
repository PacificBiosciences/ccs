//
//  MathUtils.h
//  ConsensusCoreXCode
//
//  Created by Nigel Delaney on 2/23/15.
//  Copyright (c) 2015 Pacific Biosciences. All rights reserved.
//

#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <math.h>
#include <cmath>

#define NEG_INF - INFINITY

#if YEPP
#include <ConsensusCore/Arrow/Transcendentals.hpp>
#endif

const double log_one_third = -1.098612289;

// computes log(a + b) given log(a) and log(b)
inline double logadd(double lna, double lnb)
{
    double max_val, min_val;
    if (lna >= lnb) {
        max_val = lna;
        min_val = lnb;
    }
    else {
        max_val = lnb;
        min_val = lna;
    }
    if (min_val == NEG_INF) {
        return  max_val;
    }
    min_val -= max_val;
    #if YEPP
    auto sum = yepBuiltin_Exp_64f_64f(min_val) + 1.0;
    #else
    auto sum = std::exp(min_val) + 1.0;
    #endif
    
    sum = std::log(sum) + max_val;
    return sum;
}
// Computes the log of the summed exp values for these.
inline double logsumlog(double v1, double v2, double v3, double v4)
{
    auto max_val = std::max(std::max(std::max(v1, v2), v3),v4);
    if (max_val == NEG_INF)
    {
        return max_val;
    }
    v1 -= max_val;
    v2 -= max_val;
    v3 -= max_val;
    v4 -= max_val;
    auto sum = exp(v1) + exp (v2) + exp(v3) + exp (v4);
    return log (sum) + max_val;
}
