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

#include <boost/range/as_array.hpp>
#include <boost/tuple/tuple.hpp>
#include <string>
#include <utility>
#include <vector>

#include <ConsensusCore/Mutation.hpp>
#include <ConsensusCore/Types.hpp>
#include <ConsensusCore/Utils.hpp>
#include <ConsensusCore/Arrow/MutationEnumerator.hpp>

namespace ConsensusCore {
namespace Arrow {

    namespace { // PRIVATE
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
    } // PRIVATE


    namespace detail {

        AbstractMutationEnumerator::AbstractMutationEnumerator(const std::string& tpl)
        : tpl_(tpl)
        {}

        AbstractMutationEnumerator::~AbstractMutationEnumerator() {}

    } // detail


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


    DinucleotideRepeatMutationEnumerator::DinucleotideRepeatMutationEnumerator(const std::string& tpl,
                                                                 int minDinucRepeatElements)
        : detail::AbstractMutationEnumerator(tpl)
        , minDinucRepeatElements_(minDinucRepeatElements)
    {}

    std::vector<Mutation>
    DinucleotideRepeatMutationEnumerator::Mutations() const
    {
        return Mutations(0, tpl_.length());
    }

    std::vector<Mutation>
    DinucleotideRepeatMutationEnumerator::Mutations(int beginPos, int endPos) const
    {
        std::vector<Mutation> result;

        if (minDinucRepeatElements_ <= 0)
            return result;

        //
        // Consider all dinucleotide repeats that _start_ in the window
        // and don't increment pos here, that happens later
        //
        boost::tie(beginPos, endPos) = BoundInterval(tpl_, beginPos, endPos);
        for (int pos = beginPos; pos + 1 < endPos;)
        {
            char x = tpl_[pos];
            char y = tpl_[pos + 1];
            int numElements = 1;

            for (int i = pos + 2; i + 1 < static_cast<int>(tpl_.length()); i += 2)
            {
                //
                // if we're beyond the window and have enough elements, stop
                //
                if (numElements >= minDinucRepeatElements_ && i >= endPos)
                    break;

                else if (tpl_[i] == x && tpl_[i + 1] == y)
                    numElements++;
                else break;
            }

            if (numElements >= minDinucRepeatElements_)
            {
                std::string dinuc;
                dinuc.push_back(x);
                dinuc.push_back(y);
                result.push_back(Mutation(INSERTION, pos, pos, dinuc));
                result.push_back(Mutation(DELETION, pos, pos + 2, std::string()));
            }

            //
            // Increment by the number of repeats found (if > 1),
            // less 1 so that we resume on the last base of the repeat,
            // otherwise just move ahead by 1
            //
            if (numElements > 1)
                pos += 2 * numElements - 1;
            else pos++;
        }

        return result;
    }
}
}
