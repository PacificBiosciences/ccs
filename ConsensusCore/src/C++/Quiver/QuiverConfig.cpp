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

#include <ConsensusCore/Quiver/QuiverConfig.hpp>

#include <stdexcept>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>

namespace ConsensusCore {
    QuiverConfig::QuiverConfig(const QvModelParams& qvParams,
                               int movesAvailable,
                               const BandingOptions& bandingOptions,
                               float fastScoreThreshold,
                               float addThreshold)
        : QvParams(qvParams),
          MovesAvailable(movesAvailable),
          Banding(bandingOptions),
          FastScoreThreshold(fastScoreThreshold),
          AddThreshold(addThreshold)
    {}

    QuiverConfig::QuiverConfig(const QuiverConfig& qvConfig)
        : QvParams(qvConfig.QvParams),
          MovesAvailable(qvConfig.MovesAvailable),
          Banding(qvConfig.Banding),
          FastScoreThreshold(qvConfig.FastScoreThreshold),
          AddThreshold(qvConfig.AddThreshold)
    {}


    QuiverConfigTable::QuiverConfigTable()
    {}

    bool QuiverConfigTable::InsertAs_(const std::string& name, const QuiverConfig& config)
    {
        const_iterator it;

        for (it = table.begin(); it != table.end(); it++)
            if (name.compare(it->first) == 0)
                return false;

        table.push_front(std::make_pair(name, config));

        return true;
    }

    #define FALLBACK "*"

    bool QuiverConfigTable::InsertAs(const std::string& name, const QuiverConfig& config)
        throw(InvalidInputError)

    {
        if (name.compare(FALLBACK) == 0)
            throw InvalidInputError(
                "Cannot Insert(...) a QuiverConfig with chemistry '" FALLBACK "'");
        return InsertAs_(name, config);
    }

    bool QuiverConfigTable::Insert(const QuiverConfig& config)
        throw(InvalidInputError)
    {
        const std::string& name = config.QvParams.ChemistryName;
        return InsertAs(name, config);
    }

    bool QuiverConfigTable::InsertDefault(const QuiverConfig& config)
    {
        return InsertAs_(FALLBACK, config);
    }

    int QuiverConfigTable::Size() const
    {
        return table.size();
    }

    const QuiverConfig& QuiverConfigTable::At(const std::string& name) const
        throw(InvalidInputError)
    {
        const_iterator it;

        // If we find a direct match for the chemistry, use it
        for (it = table.begin(); it != table.end(); it++)
            if (it->first.compare(name) == 0)
                return it->second;

        // Fallback
        for (it = table.begin(); it != table.end(); it++)
            if (it->first.compare(FALLBACK) == 0)
                return it->second;

        throw InvalidInputError("Chemistry not found in QuiverConfigTable");
    }


    std::vector<std::string> QuiverConfigTable::Keys() const
    {
        std::vector<std::string> keys;
        for (const_iterator it = table.begin(); it != table.end(); it++)
        {
            keys.push_back(it->first);
        }
        return keys;
    }

    QuiverConfigTable::const_iterator QuiverConfigTable::begin() const
    {
        return table.begin();
    }

    QuiverConfigTable::const_iterator QuiverConfigTable::end() const
    {
        return table.end();
    }
}
