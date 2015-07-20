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

#include <boost/range.hpp>
#include <boost/shared_array.hpp>
#include <boost/utility.hpp>
#include <string>
#include <vector>

#include <ConsensusCore/Feature.hpp>
#include <ConsensusCore/Types.hpp>

namespace ConsensusCore
{
    /// \brief An object containing observed features from a sequencing run.
    struct SequenceFeatures
    {
    public:
        explicit SequenceFeatures(const std::string& seq);
        int Length() const             { return sequence_.Length(); }
        Feature<char> Sequence() const { return sequence_; }

        /// Access to the sequence bases
        const char& operator[] (int i) const { return sequence_[i]; }
        char ElementAt(int i) const          { return (*this)[i]; }


    private:
        Feature<char> sequence_;
    };

    /// \brief A features object that contains PulseToBase QV metrics
    struct QvSequenceFeatures : public SequenceFeatures
    {
        Feature<float> SequenceAsFloat;
        Feature<float> InsQv;
        Feature<float> SubsQv;
        Feature<float> DelQv;
        Feature<float> DelTag;
        Feature<float> MergeQv;

        explicit QvSequenceFeatures(const std::string& seq);

        QvSequenceFeatures(const std::string& seq,
                           const float* insQv,
                           const float* subsQv,
                           const float* delQv,
                           const float* delTag,
                           const float* mergeQv);

        QvSequenceFeatures(const std::string& seq,
                           const Feature<float> insQv,
                           const Feature<float> subsQv,
                           const Feature<float> delQv,
                           const Feature<float> delTag,
                           const Feature<float> mergeQv);

        QvSequenceFeatures(const std::string& seq,
                           const unsigned char* insQv,
                           const unsigned char* subsQv,
                           const unsigned char* delQv,
                           const unsigned char* delTag,
                           const unsigned char* mergeQv);
    };

    /// \brief A features object for the Arrow model (contains IQV data)
    struct ArrowSequenceFeatures : public SequenceFeatures
    {
        Feature<unsigned char> InsQv;

        explicit ArrowSequenceFeatures(const std::string& seq);

        ArrowSequenceFeatures(const std::string& seq,
                              const unsigned char* insQv);

        ArrowSequenceFeatures(const std::string& seq,
                              const std::vector<unsigned char>& insQv);
    };

    /// \brief A features object that contains sequence in channel space.
    struct ChannelSequenceFeatures : public SequenceFeatures
    {
        Feature<int> Channel;

        explicit ChannelSequenceFeatures(const std::string& seq);

        ChannelSequenceFeatures(const std::string& seq, const std::vector<int>& channel);
    };
}
