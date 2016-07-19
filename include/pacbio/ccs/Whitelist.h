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

#include <map>
#include <stdexcept>
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/optional.hpp>

#include <pacbio/data/IntervalTree.h>

namespace PacBio {
namespace CCS {

class Whitelist
{
public:
    Whitelist(const std::string& spec) : all(false), globalZmws(boost::none)
    {
        // if we're all or *:*, then set all and break out
        if (spec == "all" || spec == "*:*") {
            all = true;
            return;
        }

        std::vector<std::string> mspecs;
        boost::split(mspecs, spec, boost::is_any_of(";"));

        for (const auto& mspec : mspecs) {
            // no craziness policy
            if (mspec == "all" || mspec == "*:*" || globalZmws)
                throw std::invalid_argument("invalid whitelist specification");

            std::vector<std::string> parts;
            boost::split(parts, mspec, boost::is_any_of(":"));

            // only 1 part, it's a ZMW range
            if (parts.size() == 1) {
                if (movieZmws.empty()) {
                    globalZmws = PacBio::Data::IntervalTree::FromString(parts[0]);
                    continue;
                }
            }
            // two parts, but *:_? again, it's just a ZMW range
            else if (parts.size() == 2 && parts[0] == "*") {
                if (movieZmws.empty()) {
                    globalZmws = PacBio::Data::IntervalTree::FromString(parts[1]);
                    continue;
                }
            }
            // two parts, either _:_ or _:*?
            //   either grab a range from the movie or everything, respectively
            else if (parts.size() == 2 && movieZmws.find(parts[0]) == movieZmws.end()) {
                if (parts[1] == "*")
                    movieZmws[parts[0]] = boost::none;
                else
                    movieZmws[parts[0]] = PacBio::Data::IntervalTree::FromString(parts[1]);

                continue;
            }

            // anything else is bad, including resetting any range
            throw std::invalid_argument("invalid whitelist specification");
        }
    }

    bool Contains(const std::string& movieName, int32_t holeNumber) const
    {
        if (all) return true;

        if (globalZmws) return globalZmws->Contains(holeNumber);

        auto it = movieZmws.find(movieName);
        if (it != movieZmws.end()) return !it->second || it->second->Contains(holeNumber);

        return false;
    }

private:
    bool all;
    boost::optional<PacBio::Data::IntervalTree> globalZmws;
    std::map<std::string, boost::optional<PacBio::Data::IntervalTree>> movieZmws;
};

}  // namespace CCS
}  // namespace PacBio
