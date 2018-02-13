// Copyright (c) 2017, Pacific Biosciences of California, Inc.
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

// Author: Derek Barnett

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include <pacbio/genomicconsensus/experimental/NoCallStyle.h>
#include <pacbio/genomicconsensus/experimental/ReferenceWindow.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

///
/// \brief The Consensus struct
///
struct Consensus
{
    ReferenceWindow window;
    std::string sequence;
    std::vector<uint8_t> confidence;

    ///
    /// \brief NoCallConsensus
    /// \param style
    /// \param window
    /// \param refSeq
    /// \return
    ///
    static Consensus NoCallConsensus(const NoCallStyle style, const ReferenceWindow& window,
                                     const std::string& refSeq);

    ///
    /// \brief Join
    /// \param subconsensi
    /// \return
    ///
    static Consensus Join(std::vector<Consensus> subconsensi);
};

inline bool operator<(const Consensus& lhs, const Consensus& rhs)
{
    return lhs.window < rhs.window;
}
inline bool operator==(const Consensus& lhs, const Consensus& rhs)
{
    return lhs.window == rhs.window;
}
inline bool operator!=(const Consensus& lhs, const Consensus& rhs)
{
    return lhs.window != rhs.window;
}

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
