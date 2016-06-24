// Copyright (c) 2011-2016, Pacific Biosciences of California, Inc.
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

// Authors: David Alexander, Lance Hepler

#include <algorithm>
#include <cassert>
#include <cfloat>
#include <string>
#include <vector>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <pacbio/consensus/Sequence.h>
#include <pacbio/consensus/align/AffineAlignment.h>
#include <pacbio/consensus/align/PairwiseAlignment.h>
// #include <pacbio/consensus/Utils.h>

namespace PacBio {
namespace Consensus {
namespace {

class IupacAware;
class Standard;

static float MAX4(float a, float b, float c, float d)
{
    return std::max(std::max(a, b), std::max(c, d));
}

inline bool IsIupacPartialMatch(char iupacCode, char b)
{
    assert(iupacCode != b);

    switch (iupacCode) {
        case 'R':
            return (b == 'A' || b == 'G');
        case 'Y':
            return (b == 'C' || b == 'T');
        case 'S':
            return (b == 'G' || b == 'C');
        case 'W':
            return (b == 'A' || b == 'T');
        case 'K':
            return (b == 'G' || b == 'T');
        case 'M':
            return (b == 'A' || b == 'C');
        default:
            return false;
    }
}

template <typename C>
inline float MatchScore(char t, char q, float matchScore, float mismatchScore,
                        float partialMatchScore);

template <>
inline float MatchScore<Standard>(char t, char q, float matchScore, float mismatchScore,
                                  float partialMatchScore)
{
    return (t == q ? matchScore : mismatchScore);
}

template <>
inline float MatchScore<IupacAware>(char t, char q, float matchScore, float mismatchScore,
                                    float partialMatchScore)
{
    if (t == q) {
        return matchScore;
    } else if (IsIupacPartialMatch(t, q)) {
        return partialMatchScore;
    } else if (IsIupacPartialMatch(q, t)) {
        return partialMatchScore;
    } else {
        return mismatchScore;
    }  // NOLINT
}

template <class C>
PairwiseAlignment* AlignAffineGeneric(const std::string& target, const std::string& query,
                                      AffineAlignmentParams params)
{
    // Implementation follows the textbook "two-state" affine gap model
    // description from Durbin et. al
    using boost::numeric::ublas::matrix;

    int I = query.length();
    int J = target.length();
    matrix<float> M(I + 1, J + 1);
    matrix<float> GAP(I + 1, J + 1);

    // Initialization
    M(0, 0) = 0;
    GAP(0, 0) = -FLT_MAX;
    for (int i = 1; i <= I; ++i) {
        M(i, 0) = -FLT_MAX;
        GAP(i, 0) = params.GapOpen + (i - 1) * params.GapExtend;
    }
    for (int j = 1; j <= J; ++j) {
        M(0, j) = -FLT_MAX;
        GAP(0, j) = params.GapOpen + (j - 1) * params.GapExtend;
    }

    // Main part of the recursion
    for (int i = 1; i <= I; ++i) {
        for (int j = 1; j <= J; ++j) {
            float matchScore = MatchScore<C>(target[j - 1], query[i - 1], params.MatchScore,
                                             params.MismatchScore, params.PartialMatchScore);
            M(i, j) = std::max(M(i - 1, j - 1), GAP(i - 1, j - 1)) + matchScore;
            GAP(i, j) = MAX4(M(i, j - 1) + params.GapOpen, GAP(i, j - 1) + params.GapExtend,
                             M(i - 1, j) + params.GapOpen, GAP(i - 1, j) + params.GapExtend);
        }
    }

    // Perform the traceback
    const int MATCH_MATRIX = 1;
    const int GAP_MATRIX = 2;

    std::string raQuery, raTarget;
    int i = I, j = J;
    int mat = (M(I, J) >= GAP(I, J) ? MATCH_MATRIX : GAP_MATRIX);
    int iPrev, jPrev, matPrev;
    while (i > 0 || j > 0) {
        if (mat == MATCH_MATRIX) {
            matPrev = (M(i - 1, j - 1) >= GAP(i - 1, j - 1) ? MATCH_MATRIX : GAP_MATRIX);
            iPrev = i - 1;
            jPrev = j - 1;
            raQuery.push_back(query[iPrev]);
            raTarget.push_back(target[jPrev]);
        } else {
            assert(mat == GAP_MATRIX);
            float s[4];
            s[0] = (j > 0 ? M(i, j - 1) + params.GapOpen : -FLT_MAX);
            s[1] = (j > 0 ? GAP(i, j - 1) + params.GapExtend : -FLT_MAX);
            s[2] = (i > 0 ? M(i - 1, j) + params.GapOpen : -FLT_MAX);
            s[3] = (i > 0 ? GAP(i - 1, j) + params.GapExtend : -FLT_MAX);
            int argMax = std::max_element(s, s + 4) - s;

            matPrev = ((argMax == 0 || argMax == 2) ? MATCH_MATRIX : GAP_MATRIX);
            if (argMax == 0 || argMax == 1) {
                iPrev = i;
                jPrev = j - 1;
                raQuery.push_back('-');
                raTarget.push_back(target[jPrev]);
            } else {
                iPrev = i - 1;
                jPrev = j;
                raQuery.push_back(query[iPrev]);
                raTarget.push_back('-');
            }
        }

        // Go to previous square
        i = iPrev;
        j = jPrev;
        mat = matPrev;
    }

    assert(raQuery.length() == raTarget.length());
    return new PairwiseAlignment(Reverse(raTarget), Reverse(raQuery));
}

}  // anonymous namespace

AffineAlignmentParams::AffineAlignmentParams(float matchScore, float mismatchScore, float gapOpen,
                                             float gapExtend, float partialMatchScore)
    : MatchScore(matchScore)
    , MismatchScore(mismatchScore)
    , GapOpen(gapOpen)
    , GapExtend(gapExtend)
    , PartialMatchScore(partialMatchScore)
{
}

AffineAlignmentParams DefaultAffineAlignmentParams()
{
    return AffineAlignmentParams(0, -1.0, -1.0, -0.5, 0);
}

AffineAlignmentParams IupacAwareAffineAlignmentParams()
{
    return AffineAlignmentParams(0, -1.0, -1.0, -0.5, -0.25);
}

PairwiseAlignment* AlignAffine(const std::string& target, const std::string& query,
                               AffineAlignmentParams params)
{
    return AlignAffineGeneric<Standard>(target, query, params);
}

PairwiseAlignment* AlignAffineIupac(const std::string& target, const std::string& query,
                                    AffineAlignmentParams params)
{
    return AlignAffineGeneric<IupacAware>(target, query, params);
}

}  // namespace Consensus
}  // namespace PacBio
