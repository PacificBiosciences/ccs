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

// Author: Lance Hepler, Brett Bowman

#pragma once

#include <vector>

#include <seqan/basic.h>

namespace PacBio {
namespace SparseAlignment {

template<typename TShape>
class HomopolymerHasher
{

public:  // Structors

    // Default constructor
    HomopolymerHasher(TShape& shape)
    {
        Initialize(shape);
    }

    // Move constructor
    HomopolymerHasher(HomopolymerHasher&& src) = delete;
    // Copy constructor
    HomopolymerHasher(const HomopolymerHasher& src) = delete;
    // Move assignment operator
    HomopolymerHasher& operator=(HomopolymerHasher&& rhs) = delete;
    // Copy assignment operator
    HomopolymerHasher& operator=(const HomopolymerHasher& rhs) = delete;
    // Destructor
    ~HomopolymerHasher() = default;

private:  // Internal class methods

    /// Initialize the hasher object by creating all possible
    /// homopolymer templates that fill the shape used by the
    /// index and saves them for later used by 
    ///
    /// \param  shape  The Q-gram shape used in the index
    void Initialize(const TShape& shape)
    {
        using namespace seqan;

        const char dna[4] = {'A', 'C', 'G', 'T'};

        for (size_t i = 0; i < 4; i++)
        {
            DnaString s = std::string(length(shape), dna[i]);
            hashes[i] = seqan::hash(shape, begin(s));
        }
    }

public:  // Main Function

    /// Given the hash of a Q-Gram, check whether it matches the
    /// pattern of one of the stored hashes for known homopolymers.
    /// If so return True, otherwise False.
    /// 
    /// \param  h  The Q-Gram hash of the query sequence
    inline
    bool operator()(const unsigned h) const
    {
        if (h == hashes[0] || h == hashes[1] ||
            h == hashes[2] || h == hashes[3])
            return true;

        return false;
    }

private:  // Class variables
    unsigned hashes[4];

};  // class HomopolymerHasher

}  // namespace SparseAlignment
}  // namespace PacBio
