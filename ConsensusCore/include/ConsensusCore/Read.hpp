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

#include <string>

#include <ConsensusCore/Features.hpp>
#include <ConsensusCore/Types.hpp>

namespace ConsensusCore
{
    template<typename FeaturesType>
    struct Read
    {
        FeaturesType Features;
        std::string Name;
        std::string Chemistry;

        Read(FeaturesType features,
             std::string name,
             std::string chemistry);

        Read(const Read& other);

        int Length() const;
        std::string ToString() const;

        static Read Null();
    };

    enum StrandEnum
    {
        FORWARD_STRAND = 0,
        REVERSE_STRAND = 1
    };

    template<typename FeaturesType>
    struct MappedRead : public Read<FeaturesType>
    {
        StrandEnum Strand;
        int TemplateStart;
        int TemplateEnd;
        bool PinStart;
        bool PinEnd;

        MappedRead(const Read<FeaturesType>& read,
                   StrandEnum strand,
                   int templateStart,
                   int templateEnd,
                   bool pinStart = true,
                   bool pinEnd = true);

        MappedRead(const MappedRead<FeaturesType>& other);

        std::string ToString() const;
    };

    typedef Read<QvSequenceFeatures>    QvRead;
    typedef Read<ArrowSequenceFeatures> ArrowRead;

    typedef MappedRead<QvSequenceFeatures>    MappedQvRead;
    typedef MappedRead<ArrowSequenceFeatures> MappedArrowRead;
}
