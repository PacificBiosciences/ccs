// Copyright (c) 2016, Pacific Biosciences of California, Inc.
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

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <vector>

#include <pbbam/EntireFileQuery.h>
#include <pbbam/PbiFilterQuery.h>

#include <pacbio/data/ArrayRead.h>
#include <pacbio/data/MSA.h>
#include <pacbio/io/Utility.h>
#include <pbbam/BamRecord.h>

#include <pacbio/realign/Fuse.h>

namespace PacBio {
namespace Realign {

Fuse::Fuse(const std::string& ccsInput)
{
    const auto arrayReads = FetchAlignedReads(ccsInput);
    consensusSequence_ = CreateConsensus(arrayReads);
}
Fuse::Fuse(const std::vector<Data::ArrayRead>& arrayReads)
{
    consensusSequence_ = CreateConsensus(arrayReads);
}

std::string Fuse::CreateConsensus(const std::vector<Data::ArrayRead>& arrayReads) const
{
    Data::MSA msa(arrayReads);

    auto posInsCov = CollectInsertions(msa);
    std::map<int, std::string> posIns;
    while (!posInsCov.empty())
        posIns.insert(FindInsertions(&posInsCov));

    std::string consensus;
    for (const auto& c : msa.counts) {
        if (posIns.find(c.refPos) != posIns.cend()) consensus += posIns[c.refPos];
        if (c.Coverage() > minCoverage_) {
            const auto maxBase = c.MaxBase();
            if (maxBase != '-') consensus += c.MaxBase();
        }
    }
    return consensus;
}

std::map<int, std::pair<std::string, int>> Fuse::CollectInsertions(const Data::MSA& msa) const
{
    std::map<int, std::pair<std::string, int>> posInsCov;
    for (const auto& c : msa) {
        if (!c.insertions.empty()) {
            int argmax = -1;
            std::string max;
            for (const auto& ins_count : c.insertions) {
                if (ins_count.first.size() % 3 != 0) continue;
                if (ins_count.second > argmax && ins_count.second > minInsertionCoverage_) {
                    argmax = ins_count.second;
                    max = ins_count.first;
                }
            }
            if (argmax != -1) posInsCov[c.refPos] = std::make_pair(max, argmax);
        }
    }
    return posInsCov;
}

std::pair<int, std::string> Fuse::FindInsertions(
    std::map<int, std::pair<std::string, int>>* posInsCov, int windowSize) const
{
    int argMax = -1;
    int maxCov = -1;
    std::string ins;
    for (const auto& kv : *posInsCov) {
        if (kv.second.second > maxCov) {
            maxCov = kv.second.second;
            argMax = kv.first;
            ins = kv.second.first;
        }
    }

    for (int i = std::max(0, argMax - windowSize); i < argMax + windowSize; ++i) {
        auto it = posInsCov->find(i);
        if (it != posInsCov->cend()) posInsCov->erase(it);
    }
    return std::make_pair(argMax, ins);
}

std::vector<Data::ArrayRead> Fuse::FetchAlignedReads(const std::string& ccsInput) const
{
    BAM::DataSet ds(ccsInput);
    const auto filter = BAM::PbiFilter::FromDataSet(ds);
    std::unique_ptr<BAM::internal::IQuery> query(nullptr);
    if (filter.IsEmpty())
        query.reset(new BAM::EntireFileQuery(ds));
    else
        query.reset(new BAM::PbiFilterQuery(filter, ds));

    std::vector<Data::ArrayRead> reads;
    int idx = 0;
    for (const auto& read : *query) {
        reads.emplace_back(Data::BAMArrayRead(read, idx++));
    }

    return reads;
}
}
}  // ::PacBio::Realign