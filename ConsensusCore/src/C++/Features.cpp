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

#include <ConsensusCore/Features.hpp>

#include <algorithm>
#include <string>
#include <vector>

#include <ConsensusCore/Feature.hpp>
#include <ConsensusCore/Types.hpp>
#include <ConsensusCore/Utils.hpp>


ConsensusCore::SequenceFeatures::SequenceFeatures(const std::string& seq)
    : sequence_(seq.c_str(), seq.length())
{}

namespace
{
    void CheckTagFeature(ConsensusCore::Feature<float> feature)
    {
        foreach (const float& tag, feature)
        {
            if (!(tag == 'A' ||
                  tag == 'C' ||
                  tag == 'G' ||
                  tag == 'T' ||
                  tag == 'N' ||
                  tag ==  0))
            {
                throw ConsensusCore::InternalError(
                    "Invalid DelTag field---must be ASCII encoded.");
            }
        }
    }
}

namespace ConsensusCore
{
    QvSequenceFeatures::QvSequenceFeatures(const std::string& seq)
        : SequenceFeatures(seq),
          SequenceAsFloat(Length()),
          InsQv (Length()),
          SubsQv(Length()),
          DelQv (Length()),
          DelTag(Length()),
          MergeQv(Length())
    {
        for (int i = 0; i < Length(); i++)
        {
            SequenceAsFloat[i] = static_cast<float>(seq[i]);
        }
    }

    QvSequenceFeatures::QvSequenceFeatures(const std::string& seq,
                                           const float* insQv,
                                           const float* subsQv,
                                           const float* delQv,
                                           const float* delTag,
                                           const float* mergeQv)
        : SequenceFeatures(seq),
          SequenceAsFloat(Length()),
          InsQv  (insQv,   Length()),
          SubsQv (subsQv,  Length()),
          DelQv  (delQv,   Length()),
          DelTag (delTag,  Length()),
          MergeQv(mergeQv, Length())
    {
        for (int i = 0; i < Length(); i++)
        {
            SequenceAsFloat[i] = static_cast<float>(seq[i]);
        }
        CheckTagFeature(DelTag);
    }


    QvSequenceFeatures::QvSequenceFeatures(const std::string& seq,
                                           const unsigned char* insQv,
                                           const unsigned char* subsQv,
                                           const unsigned char* delQv,
                                           const unsigned char* delTag,
                                           const unsigned char* mergeQv)
        : SequenceFeatures(seq),
          SequenceAsFloat(Length()),
          InsQv  (reinterpret_cast<const float*>(insQv),   Length()),
          SubsQv (reinterpret_cast<const float*>(subsQv),  Length()),
          DelQv  (reinterpret_cast<const float*>(delQv),   Length()),
          DelTag (reinterpret_cast<const float*>(delTag),  Length()),
          MergeQv(reinterpret_cast<const float*>(mergeQv), Length())
    {
        for (int i = 0; i < Length(); i++)
        {
            SequenceAsFloat[i] = static_cast<float>(seq[i]);
        }
        CheckTagFeature(DelTag);
    }


    QvSequenceFeatures::QvSequenceFeatures(const std::string& seq,
                                           const Feature<float> insQv,
                                           const Feature<float> subsQv,
                                           const Feature<float> delQv,
                                           const Feature<float> delTag,
                                           const Feature<float> mergeQv)
        : SequenceFeatures(seq),
          SequenceAsFloat(Length()),
          InsQv (insQv),
          SubsQv(subsQv),
          DelQv (delQv),
          DelTag(delTag),
          MergeQv(mergeQv)
    {
        for (int i = 0; i < Length(); i++)
        {
            SequenceAsFloat[i] = static_cast<float>(seq[i]);
        }
        CheckTagFeature(DelTag);
    }


    ArrowSequenceFeatures::ArrowSequenceFeatures(const std::string& seq)
        : SequenceFeatures(seq)
        , InsQv(Length())
    {}

    ArrowSequenceFeatures::ArrowSequenceFeatures(const std::string& seq,
                                                 const unsigned char* const insQv)
        : SequenceFeatures(seq)
        , InsQv(insQv, Length())
    {}

    ArrowSequenceFeatures::ArrowSequenceFeatures(const std::string& seq,
                                                 const std::vector<unsigned char>& insQv)
        : SequenceFeatures(seq)
        , InsQv(&(insQv[0]), Length())
    {}


    ChannelSequenceFeatures::ChannelSequenceFeatures(const std::string& seq)
        : SequenceFeatures(seq),
          Channel(Length())
    {}

    ChannelSequenceFeatures::ChannelSequenceFeatures(const std::string& seq,
                                                     const std::vector<int>& channel)
        : SequenceFeatures(seq),
          Channel(&(channel[0]), Length())
    {}
}
