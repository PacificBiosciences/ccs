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

// Author: Lance Hepler

#pragma once

#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/sequence.h>

namespace PacBio {
namespace SparseAlignment {

/// A simple struct for represting a complete specialization for a SeqAn index
/// such that template-ifying code around different indices can be done around
/// a single, simple parameter.
///
/// TODO (lhepler) : investigate default values other than 10
///
/// \param  size_t  The length of the K-mer or Q-Gram used by the index.
///                 By default 10.
/// \param  TShape  The SeqAn shape specialization.  By default Ungapped.
/// \param  TIndex  The SeqAn index specialization.  By default Q-Gram.
template<size_t   TSize  = 10,
         typename TShape = seqan::UngappedShape<TSize>,
         typename TIndex = seqan::IndexQGram<TShape>>
struct FindSeedsConfig
{
    typedef TIndex    IndexType;
    typedef TShape    ShapeType;
    static const size_t Size = TSize;
};

} // SparseAlignment
} // PacBio
