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

// Author: Armin Töpfer

#pragma once

#include <pacbio/data/FisherResult.h>
#include <pacbio/statistics/Fisher.h>
#include <algorithm>

namespace PacBio {
namespace Statistics {
class Tests
{
public:
    /// Compute Fisher's exact test for CCS substitutions and deletions
    static std::map<std::string, double> FisherCCS(const std::array<int, 5>& observed, const std::map<std::string, int> insertions, const double threshold)
    {
        int argMax = 0;
        double sum = 0;
        const auto pml = CalculatePml(observed, &argMax, &sum);

        std::map<std::string, double> results;
        for (const auto& kv : insertions)
        {
            const double p = Fisher::fisher_exact_tiss(kv.second + 1, sum, 0.0084 / 4.0 * sum, sum);
            if (p < threshold)
                results.insert({kv.first, p});
        }

        return results;
    }

    /// Compute Fisher's exact test for CCS substitutions and deletions
    static Data::FisherResult FisherCCS(const std::array<int, 5>& observed, const double threshold)
    {
        int argMax = 0;
        double sum = 0;
        const auto pml = CalculatePml(observed, &argMax, &sum);
        const auto pMatch = CalculatePriors(argMax);

        auto fisherCCS = [&observed, &pMatch, &pml, &sum](int i) {
            return Fisher::fisher_exact_tiss((pml[i] * sum), (sum), (pMatch[i] * sum), (sum));
        };

        Data::FisherResult fr;
        fr.pValues = {{fisherCCS(0), fisherCCS(1), fisherCCS(2), fisherCCS(3), fisherCCS(4)}};
        for (int i = 0; i < 5; ++i) {
            if (fr.pValues.at(i) < threshold && observed.at(i) > 1) {
                if (i != argMax) fr.hit = true;
                fr.mask[i] = 1;
            }
            // else if (i == argMax)
            //     fr.mask[i] = 1;
            else
                fr.mask[i] = 0;
        }
        fr.argMax = argMax;
        return fr;
    }

private:
    static std::array<double, 5> CalculatePml(const std::array<int, 5>& observed, int* argMax,
                                              double* sum)
    {
        std::array<double, 5> pml;
        std::copy_n(std::cbegin(observed), 5, std::begin(pml));

        // +1 each entry
        std::for_each(std::begin(pml), std::end(pml), [](double& p) { ++p; });

        *argMax =
            std::distance(std::cbegin(pml), std::max_element(std::cbegin(pml), std::cend(pml)));
        *sum = std::accumulate(std::cbegin(pml), std::cend(pml), 0.0);

        // normalize
        std::for_each(std::begin(pml), std::end(pml), [&sum](double& p) { p /= (*sum); });
        return pml;
    }

    static std::array<double, 5> CalculatePriors(const int argMax)
    {
        assert(argMax >= 0 && argMax <= 5);

        std::array<double, 5> pMatch{{0.0005, 0.0005, 0.0005, 0.0005, 0.0029}};
        pMatch[argMax] = 0.9872;
        double pMatchSum = pMatch[0] + pMatch[1] + pMatch[2] + pMatch[3] + pMatch[4];
        for (auto& p : pMatch)
            p /= pMatchSum;

        return pMatch;
    }
};
}
}  // ::PacBio::Statistics