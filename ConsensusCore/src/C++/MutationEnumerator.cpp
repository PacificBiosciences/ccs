// Copyright (c) 2011-2014, Pacific Biosciences of California, Inc.
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

#include <ConsensusCore/MutationEnumerator.hpp>
#include <ConsensusCore/Types.hpp>
#include <ConsensusCore/Utils.hpp>
#include <ConsensusCore/Mutation.hpp>

#include <boost/range/as_array.hpp>
#include <boost/tuple/tuple.hpp>
#include <cstring>
#include <string>
#include <utility>
#include <vector>

namespace ConsensusCore
{
    namespace {  // PRIVATE
        const char BASES[] = { 'A', 'C', 'G', 'T' };

        inline
        int BoundPosition(const std::string& tpl, int pos)
        {
            return pos < 0 ? 0 : (pos > static_cast<int>(tpl.length()) ? tpl.length() : pos);
        }

        inline
        std::pair<int, int> BoundInterval(const std::string& tpl, int beginPos, int endPos)
        {
            return std::make_pair(BoundPosition(tpl, beginPos), BoundPosition(tpl, endPos));
        }
    }  // PRIVATE


    namespace detail {

        AbstractMutationEnumerator::AbstractMutationEnumerator(const std::string& tpl)
        : tpl_(tpl)
        {}

        AbstractMutationEnumerator::~AbstractMutationEnumerator() {}
    }  // detail


    AllSingleBaseMutationEnumerator::AllSingleBaseMutationEnumerator(const std::string& tpl)
        : detail::AbstractMutationEnumerator(tpl)
    {}

    std::vector<Mutation>
    AllSingleBaseMutationEnumerator::Mutations() const
    {
        return Mutations(0, tpl_.length());
    }

    std::vector<Mutation>
    AllSingleBaseMutationEnumerator::Mutations(int beginPos, int endPos) const
    {
        std::vector<Mutation> result;
        boost::tie(beginPos, endPos) = BoundInterval(tpl_, beginPos, endPos);
        for (int pos = beginPos; pos < endPos; pos++)
        {
            foreach (char base, boost::as_array(BASES)) {
                if (base != tpl_[pos]) {
                    result.push_back(Mutation(SUBSTITUTION, pos, base));
                }
            }
            foreach (char base, boost::as_array(BASES)) {
                result.push_back(Mutation(INSERTION, pos, base));
            }
            result.push_back(Mutation(DELETION, pos, '-'));
        }
        return result;
    }


    UniqueSingleBaseMutationEnumerator::UniqueSingleBaseMutationEnumerator(const std::string& tpl)
        : detail::AbstractMutationEnumerator(tpl)
    {}

    std::vector<Mutation>
    UniqueSingleBaseMutationEnumerator::Mutations() const
    {
        return Mutations(0, tpl_.length());
    }

    std::vector<Mutation>
    UniqueSingleBaseMutationEnumerator::Mutations(int beginPos, int endPos) const
    {
        std::vector<Mutation> result;
        boost::tie(beginPos, endPos) = BoundInterval(tpl_, beginPos, endPos);
        for (int pos = beginPos; pos < endPos; pos++)
        {
            char prevTplBase = pos > 0 ? tpl_[pos-1] : '-';
            foreach (char base, boost::as_array(BASES)) {
                if (base != tpl_[pos]) {
                    result.push_back(Mutation(SUBSTITUTION, pos, base));
                }
            }
            // Insertions only allowed at the beginning of homopolymers
            foreach (char base, boost::as_array(BASES)) {
                if (base != prevTplBase) {
                    result.push_back(Mutation(INSERTION, pos, base));
                }
            }
            // Deletions only allowed at the beginning of homopolymers
            if (tpl_[pos] != prevTplBase) {
                result.push_back(Mutation(DELETION, pos, '-'));
            }
        }
        return result;
    }


    RepeatMutationEnumerator::RepeatMutationEnumerator
    (const std::string& tpl, int repeatLength, int minRepeatElements)
        : detail::AbstractMutationEnumerator(tpl)
        , repeatLength_(repeatLength)
        , minRepeatElements_(minRepeatElements)
    {}

    std::vector<Mutation>
    RepeatMutationEnumerator::Mutations() const
    {
        return Mutations(0, tpl_.length());
    }

    std::vector<Mutation>
    RepeatMutationEnumerator::Mutations(int beginPos, int endPos) const
    {
        std::vector<Mutation> result;

        if (minRepeatElements_ <= 0 || repeatLength_ > 31)
            return result;

        // buffer to store repeat characters
        char buf[32] = {0};

        //
        // Consider all dinucleotide repeats that _start_ in the window
        // and don't increment pos here, that happens later
        //
        boost::tie(beginPos, endPos) = BoundInterval(tpl_, beginPos, endPos);
        for (int pos = beginPos; pos + repeatLength_ <= endPos;)
        {
            strncpy(buf, tpl_.c_str() + pos, repeatLength_);

            int numElements = 1;

            for (int i = pos + repeatLength_;
                 static_cast<size_t>(i + repeatLength_) <= tpl_.length();
                 i += repeatLength_)
            {
                //
                // if we're beyond the window and have enough elements, stop
                //
                if (numElements >= minRepeatElements_ && i >= endPos)
                    break;

                else if (strncmp(buf, tpl_.c_str() + i, repeatLength_) == 0)
                    numElements++;
                else
                    break;
            }

            if (numElements >= minRepeatElements_)
            {
                std::string repeat(buf);
                result.push_back(Mutation(INSERTION, pos, pos, repeat));
                result.push_back(Mutation(DELETION, pos, pos + repeatLength_, std::string()));
            }

            //
            // Increment by the number of repeats found (if > 1),
            // less 1 so that we resume on the second base of the last repeat,
            // otherwise just move ahead by 1
            //
            if (numElements > 1)
                pos += repeatLength_ * (numElements - 1) + 1;
            else
                pos++;
        }

        return result;
    }


    DinucleotideRepeatMutationEnumerator::DinucleotideRepeatMutationEnumerator
    (const std::string& tpl, int minDinucRepeatElements)
        : RepeatMutationEnumerator(tpl, 2, minDinucRepeatElements)
    {}
}
