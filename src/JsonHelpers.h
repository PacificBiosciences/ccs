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

// Author: Lance Hepler

#pragma once

#include <boost/property_tree/ptree.hpp>

namespace PacBio {
namespace Consensus {
template <size_t I>
void ReadMatrix(double (&mat)[I], const boost::property_tree::ptree& pt)
{
    if (pt.size() != I) throw std::invalid_argument("bad size (1D)");
    size_t i = 0;
    for (const auto& item : pt) {
        mat[i] = item.second.get_value<double>();
        ++i;
    }
}

template <size_t I, size_t J>
void ReadMatrix(double (&mat)[I][J], const boost::property_tree::ptree& pt)
{
    if (pt.size() != I) throw std::invalid_argument("bad size (2D)");
    size_t i = 0;
    for (const auto& item : pt) {
        ReadMatrix<J>(mat[i], item.second);
        ++i;
    }
}

template <size_t I, size_t J, size_t K>
void ReadMatrix(double (&mat)[I][J][K], const boost::property_tree::ptree& pt)
{
    if (pt.size() != I) throw std::invalid_argument("bad size (3D)");
    size_t i = 0;
    for (const auto& item : pt) {
        ReadMatrix<J, K>(mat[i], item.second);
        ++i;
    }
}
}
}
