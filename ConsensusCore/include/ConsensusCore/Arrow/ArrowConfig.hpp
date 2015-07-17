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

#include <list>
#include <string>
#include <utility>
#include <vector>
#include <math.h>

#include <ConsensusCore/Utils.hpp>
#include <ConsensusCore/Arrow/MathUtils.hpp>
#include <ConsensusCore/Arrow/ContextParameters.hpp>


/* Hard coded mismatch probability for now.
   This is derived as the mean in PlotBinnedTraining.R */
#define MISMATCH_PROBABILITY 0.002671256

namespace ConsensusCore {
namespace Arrow {

    // private anonymous parameters
    namespace
    {
        const double MATCH_IQV_PMF[]  = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
        const double INSERT_IQV_PMF[] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
    }

    /// \brief The banding optimizations to be used by a recursor
    struct BandingOptions
    {
        double ScoreDiff;

        BandingOptions(double scoreDiff)
            : ScoreDiff(scoreDiff)
        {
            if(scoreDiff < 0) {
                throw InvalidInputError("ScoreDiff must be positive!");
            }
        }

    };


    /// \brief A parameter vector for analysis using the QV model
    struct ModelParams
    {
        double MatchIqvPmf[20];
        double InsertIqvPmf[20];
        double PrMiscall;
        double PrNotMiscall;
        double PrThirdOfMiscall;
        //
        // Constructor for single merge rate and merge rate slope
        //
        ModelParams(const double matchIqvPmf[]  = MATCH_IQV_PMF,
                    const double insertIqvPmf[] = INSERT_IQV_PMF,
                    double mismatch             = MISMATCH_PROBABILITY)
            : PrMiscall(mismatch)
            , PrNotMiscall(1.0 - mismatch)
            , PrThirdOfMiscall(mismatch / 3.0)
        {
            memcpy(MatchIqvPmf,  matchIqvPmf,  sizeof(MatchIqvPmf));
            memcpy(InsertIqvPmf, insertIqvPmf, sizeof(InsertIqvPmf));
        }
        
        // Copy constructor
        ModelParams(const ModelParams& src) = default;
        
        // Move constructor
        ModelParams(ModelParams&& src) = default;
        
        // Copy Assignment operator
        ModelParams& operator=(const ModelParams& rhs) = default;
    };


    class ArrowConfig
    {
        public:
            ModelParams MdlParams;
            ContextParameters CtxParams;
            BandingOptions Banding;
            double FastScoreThreshold;
            double AddThreshold;

            ArrowConfig(const ContextParameters& ctxParams,
                        const BandingOptions& bandingOptions,
                        double fastScoreThreshold = -12.5,
                        double addThreshold = 1.0f);

            ArrowConfig(const ArrowConfig& qvConfig);
        
            // Assuming compiler generated destructor is sufficient.
    };



    class ArrowConfigTable
    {
    private:
        typedef std::pair<const std::string, const ArrowConfig> ArrowConfigTableEntry;
        std::list<ArrowConfigTableEntry> table;

    public:
        typedef std::list<ArrowConfigTableEntry>::const_iterator const_iterator;

        ArrowConfigTable();

        bool Insert(const std::string& name, const ArrowConfig& config);
        int Size() const;

        const ArrowConfig& At(const std::string& name) const throw(InvalidInputError);

        std::vector<std::string> Keys() const;

#ifndef SWIG
        const_iterator begin() const;
        const_iterator end() const;
#endif
    };
}
}
