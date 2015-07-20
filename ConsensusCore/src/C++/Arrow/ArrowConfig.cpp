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

#include <stdexcept>
#include <string>
#include <vector>

#include <ConsensusCore/Arrow/ArrowConfig.hpp>

namespace ConsensusCore {
namespace Arrow {

    ArrowConfig::ArrowConfig(const ContextParameters& ctxParams,
                             const BandingOptions& bandingOptions,
                             double fastScoreThreshold,
                             double addThreshold)
        : MdlParams()
        , CtxParams(ctxParams)
        , Banding(bandingOptions)
        , FastScoreThreshold(fastScoreThreshold)
        , AddThreshold(addThreshold)
    {

    }

    ArrowConfig::ArrowConfig(const ArrowConfig& qvConfig)
        : MdlParams(qvConfig.MdlParams),
          CtxParams(qvConfig.CtxParams),
          Banding(qvConfig.Banding),
          FastScoreThreshold(qvConfig.FastScoreThreshold),
          AddThreshold(qvConfig.AddThreshold)
    {}


    ArrowConfigTable::ArrowConfigTable()
    {}

    bool ArrowConfigTable::Insert(const std::string& name, const ArrowConfig& config)
    {
        const_iterator it;

        for (it = table.begin(); it != table.end(); it++)
            if (name.compare(it->first) == 0)
                return false;

        table.push_front(std::make_pair(name, config));

        return true;
    }

    int ArrowConfigTable::Size() const
    {
        return (int)table.size();
    }

    const ArrowConfig& ArrowConfigTable::At(const std::string& name) const
        throw(InvalidInputError)
    {
        const_iterator it;

        // If we find a direct match for the chemistry, use it
        for (it = table.begin(); it != table.end(); it++)
            if (it->first.compare(name) == 0)
                return it->second;

        // Fallback is "*"
        for (it = table.begin(); it != table.end(); it++)
            if (it->first.compare("*") == 0)
                return it->second;

        throw InvalidInputError("Chemistry not found in ArrowConfigTable");
    }


    std::vector<std::string> ArrowConfigTable::Keys() const
    {
        std::vector<std::string> keys;
        for (const_iterator it = table.begin(); it != table.end(); it++)
        {
            keys.push_back(it->first);
        }
        return keys;
    }

    ArrowConfigTable::const_iterator ArrowConfigTable::begin() const
    {
        return table.begin();
    }

    ArrowConfigTable::const_iterator ArrowConfigTable::end() const
    {
        return table.end();
    }
}
}
