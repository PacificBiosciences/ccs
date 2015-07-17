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


#include <xmmintrin.h>
#include <pmmintrin.h>

#include <algorithm>
#include <cassert>
#include <cfloat>
#include <climits>
#include <cmath>
#include <iostream>
#include <string>
#include <utility>

#include <ConsensusCore/Quiver/detail/SseMath.hpp>
#include <ConsensusCore/Quiver/QuiverConfig.hpp>
#include <ConsensusCore/Features.hpp>
#include <ConsensusCore/Types.hpp>
#include <ConsensusCore/Utils.hpp>
#include <ConsensusCore/Read.hpp>

#ifndef SWIG
using std::min;
using std::max;
#endif  // SWIG

#define NEG_INF -FLT_MAX

namespace ConsensusCore
{
    //
    // Utility functions
    //
    static inline int encodeTplBase(char base)
    {
        switch (base) {
            case 'A': return 0;
            case 'C': return 1;
            case 'G': return 2;
            case 'T': return 3;
            case 'M': return 4;  // For testing
            case 'N': return 5;  // For testing
            default:  ShouldNotReachHere();
        }
    }

    //
    // Evaluator classes
    //

    /// \brief An Evaluator that can compute move scores using a QvSequenceFeatures
    class QvEvaluator
    {
    public:
        typedef QvSequenceFeatures FeaturesType;
        typedef QvModelParams      ParamsType;

    public:
        QvEvaluator(const QvRead& read,
                    const std::string& tpl,
                    const QvModelParams& params,
                    bool pinStart = true,
                    bool pinEnd = true)
            : read_(read),
              params_(params),
              tpl_(tpl),
              pinStart_(pinStart),
              pinEnd_(pinEnd)
        {}

        ~QvEvaluator()
        {}

        std::string ReadName() const
        {
            return read_.Name;
        }

        std::string Basecalls() const
        {
            return Features().Sequence();
        }

        std::string Template() const
        {
            return tpl_;
        }

        void Template(std::string tpl)
        {
            tpl_ = tpl;
        }


        int ReadLength() const
        {
            return Features().Length();
        }

        int TemplateLength() const
        {
            return tpl_.length();
        }

        bool PinEnd() const
        {
            return pinEnd_;
        }

        bool PinStart() const
        {
            return pinStart_;
        }

        bool IsMatch(int i, int j) const
        {
            assert(0 <= i && i < ReadLength());
            assert (0 <= j && j < TemplateLength());
            return (Features()[i] == tpl_[j]);
        }

        float Inc(int i, int j) const
        {
            assert(0 <= j && j < TemplateLength() &&
                   0 <= i && i < ReadLength() );
            return (IsMatch(i, j)) ?
                    params_.Match :
                    params_.Mismatch + params_.MismatchS * Features().SubsQv[i];
        }

        float Del(int i, int j) const
        {
            assert(0 <= j && j < TemplateLength() &&
                   0 <= i && i <= ReadLength() );
            if ( (!PinStart() && i == 0) || (!PinEnd() && i == ReadLength()) )
            {
                return 0.0f;
            }
            else
            {
                float tplBase = tpl_[j];
                return (i < ReadLength() && tplBase == Features().DelTag[i]) ?
                        params_.DeletionWithTag + params_.DeletionWithTagS * Features().DelQv[i] :
                        params_.DeletionN;
            }
        }

        float Extra(int i, int j) const
        {
            assert(0 <= j && j <= TemplateLength() &&
                   0 <= i && i < ReadLength() );
            return (j < TemplateLength() && IsMatch(i, j)) ?
                    params_.Branch + params_.BranchS * Features().InsQv[i] :
                    params_.Nce + params_.NceS * Features().InsQv[i];
        }

        float Merge(int i, int j) const
        {
            assert(0 <= j && j < TemplateLength() - 1 &&
                   0 <= i && i < ReadLength() );
            if (!(Features()[i] == tpl_[j] && Features()[i] == tpl_[j + 1]) )
            {
                return -FLT_MAX;
            }
            else
            {   int tplBase = encodeTplBase(tpl_[j]);
                return params_.Merge[tplBase] + params_.MergeS[tplBase] * Features().MergeQv[i];
            }
        }

        //
        // SSE
        //

        __m128 Inc4(int i, int j) const
        {
            assert (0 <= i && i <= ReadLength() - 4);
            assert (0 <= j && j < TemplateLength());
            float tplBase = tpl_[j];
            __m128 match = _mm_set_ps1(params_.Match);
            __m128 mismatch = AFFINE4(params_.Mismatch, params_.MismatchS, &Features().SubsQv[i]);
            // Mask to see it the base is equal to the template
            __m128 mask = _mm_cmpeq_ps(_mm_loadu_ps(&Features().SequenceAsFloat[i]),
                                       _mm_set_ps1(tplBase));
            return MUX4(mask, match, mismatch);
        }

        __m128 Del4(int i, int j) const
        {
            assert (0 <= i && i <= ReadLength());
            assert (0 <= j && j < TemplateLength());
            if (i != 0 && i + 3 != ReadLength())
            {
                float tplBase = tpl_[j];
                __m128 delWTag = AFFINE4(params_.DeletionWithTag,
                                         params_.DeletionWithTagS,
                                         &Features().DelQv[i]);
                __m128 delNoTag = _mm_set_ps1(params_.DeletionN);
                __m128 mask = _mm_cmpeq_ps(_mm_loadu_ps(&Features().DelTag[i]),
                                           _mm_set_ps1(tplBase));
                return MUX4(mask, delWTag, delNoTag);
            }
            else
            {
                // Have to do PinStart/PinEnd logic, and weird
                // logic for last row.  Punt.
                __m128 res = _mm_setr_ps(Del(i + 0, j),
                                         Del(i + 1, j),
                                         Del(i + 2, j),
                                         Del(i + 3, j));
                return res;
            }
        }

        __m128 Extra4(int i, int j) const
        {
            assert (0 <= i && i <= ReadLength() - 4);
            assert (0 <= j && j <= TemplateLength());
            if (i != 0 && i + 3 != ReadLength())
            {
                float tplBase = tpl_[j];
                __m128 branch = AFFINE4(params_.Branch, params_.BranchS, &Features().InsQv[i]);
                __m128 nce    = AFFINE4(params_.Nce,    params_.NceS,    &Features().InsQv[i]);

                __m128 mask = _mm_cmpeq_ps(_mm_loadu_ps(&Features().SequenceAsFloat[i]),
                                           _mm_set_ps1(tplBase));
                return MUX4(mask, branch, nce);
            }
            else
            {
                __m128 res = _mm_setr_ps(Extra(i + 0, j),
                                         Extra(i + 1, j),
                                         Extra(i + 2, j),
                                         Extra(i + 3, j));
                return res;
            }
        }


        __m128 Merge4(int i, int j) const
        {
            assert(0 <= i && i <= ReadLength() - 4);
            assert(0 <= j && j < TemplateLength() - 1);

            float tplBase     = tpl_[j];
            float tplBaseNext = tpl_[j + 1];
            int tplBase_ = encodeTplBase(tpl_[j]);

            __m128 merge =  AFFINE4(params_.Merge[tplBase_],
                                    params_.MergeS[tplBase_],
                                    &Features().MergeQv[i]);
            __m128 noMerge = _mm_set_ps1(-FLT_MAX);

            if (tplBase == tplBaseNext)
            {
                __m128 mask = _mm_cmpeq_ps(_mm_loadu_ps(&Features().SequenceAsFloat[i]),
                                           _mm_set_ps1(tplBase));
                return MUX4(mask, merge, noMerge);
            }
            else
            {
                return noMerge;
            }
        }

    protected:
        inline const QvSequenceFeatures& Features() const
        {
            return read_.Features;
        }


    protected:
        QvRead read_;
        QvModelParams params_;
        std::string tpl_;
        bool pinStart_;
        bool pinEnd_;
    };
}
