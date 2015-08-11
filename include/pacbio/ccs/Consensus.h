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

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include <OptionParser.h>

#include <ConsensusCore/Arrow/MultiReadMutationScorer.hpp>
#include <ConsensusCore/Poa/PoaConsensus.hpp>
#include <ConsensusCore/Consensus.hpp>

#include <pbbam/Accuracy.h>
#include <pbbam/QualityValues.h>
#include <pbbam/LocalContextFlags.h>

#include <pacbio/ccs/Logging.h>
#include <pacbio/ccs/SparsePoa.h>
#include <pacbio/ccs/Timer.h>


namespace PacBio {
namespace CCS {


using SNR = ConsensusCore::Arrow::SNR;
using QualityValues = PacBio::BAM::QualityValues;
using LocalContextFlags = PacBio::BAM::LocalContextFlags;
using Accuracy = PacBio::BAM::Accuracy;


namespace OptionNames {
    // constexpr auto MaxPoaCoverage       = "maxPoaCoverage";
    constexpr auto MinLength            = "minLength";
    constexpr auto MinPasses            = "minPasses";
    constexpr auto MinPredictedAccuracy = "minPredictedAccuracy";
    constexpr auto MinZScore            = "minZScore";
    constexpr auto MaxDropFraction      = "maxDropFraction";
    // constexpr auto Directional          = "directional";
} // namespace OptionNames


struct ConsensusSettings
{
    size_t MaxPoaCoverage;
    size_t MinLength;
    size_t MinPasses;
    double MinPredictedAccuracy;
    double MinZScore;
    double MaxDropFraction;
    bool   Directional;

    ConsensusSettings(const optparse::Values& options);

    static
    void AddOptions(optparse::OptionParser * const parser)
    {
        const std::string em = "--";
        // TODO(lhepler) implement alignment to POA and directional support
        // parser->add_option(em + OptionNames::MaxPoaCoverage).type("int").set_default(1024).help("Maximum number of subreads to use when building POA. Default = %default");
        parser->add_option(em + OptionNames::MinLength).type("int").set_default(10).help("Minimum length of subreads to use for generating CCS. Default = %default");
        parser->add_option(em + OptionNames::MinPasses).type("int").set_default(3).help("Minimum number of subreads required to generate CCS. Default = %default");
        parser->add_option(em + OptionNames::MinPredictedAccuracy).type("float").set_default(0.90).help("Minimum predicted accuracy in percent. Default = %default");
        parser->add_option(em + OptionNames::MinZScore).type("float").set_default(-5.0).help("Minimum z-score to use a subread. NaN disables this filter. Default = %default");
        parser->add_option(em + OptionNames::MaxDropFraction).type("float").set_default(0.33).help("Maximum fraction of subreads that can be dropped before giving up. Default = %default");
        // parser->add_option(em + OptionNames::Directional).action("store_true").set_default("0").help("Generate a consensus for each strand. Default = false");
    }
};


template<typename TId>
struct ReadType
{
    TId Id;
    std::string Seq;
    LocalContextFlags Flags;
    Accuracy ReadAccuracy;
    // TODO (move SNR here, eventually)
    // SNR SignalToNoise;
};


template<typename TId, typename TRead>
struct ChunkType
{
    TId Id;
    std::vector<TRead> Reads;
    SNR SignalToNoise;
};


template<typename TId>
struct ConsensusType
{
    TId Id;
    std::string Sequence;
    std::string Qualities;
    size_t NumPasses;
    double PredictedAccuracy;
    double GlobalZScore;
    double AvgZScore;
    std::vector<double> ZScores;
    std::vector<int32_t> StatusCounts;
    size_t MutationsTested;
    size_t MutationsApplied;
    SNR SignalToNoise;
    float ElapsedMilliseconds;
};


template<typename TConsensus>
class ResultType : public std::vector<TConsensus>
{
public:
    size_t Success;
    size_t PoorSNR;
    size_t NoSubreads;
    size_t TooShort;
    size_t TooFewPasses;
    size_t TooManyUnusable;
    size_t NonConvergent;
    size_t PoorQuality;
    size_t Other;

    ResultType()
        : Success{0}
        , PoorSNR{0}
        , NoSubreads{0}
        , TooShort{0}
        , TooFewPasses{0}
        , TooManyUnusable{0}
        , NonConvergent{0}
        , PoorQuality{0}
        , Other{0}
    { }

    ResultType<TConsensus>&
    operator+=(const ResultType<TConsensus>& other)
    {
        Success         += other.Success;
        PoorSNR         += other.PoorSNR;
        NoSubreads      += other.NoSubreads;
        TooShort        += other.TooShort;
        TooManyUnusable += other.TooManyUnusable;
        TooFewPasses    += other.TooFewPasses;
        NonConvergent   += other.NonConvergent;
        PoorQuality     += other.PoorQuality;
        Other           += other.Other;
        return *this;
    }

    size_t Total() const
    {
        return
            ( Success
            + PoorSNR
            + NoSubreads
            + TooShort
            + TooManyUnusable
            + TooFewPasses
            + NonConvergent
            + PoorQuality
            + Other );
    }
};


namespace { // anonymous

template<typename T>
float Median(std::vector<T>* vs)
{
    const size_t n = vs->size();
    std::sort(vs->begin(), vs->end());
    if (n % 2 == 1)
        return static_cast<float>(vs->at(n / 2));
    return 0.5 * (vs->at(n / 2 - 1) + vs->at(n / 2));
}

template<typename TRead>
std::vector<const TRead*> FilterReads(const std::vector<TRead>& reads,
                                      const size_t minLength)
{
    std::vector<const TRead*> results;

    if (reads.empty())
        return results;

    std::vector<size_t> lengths;

    for (const auto& read : reads)
    {
        lengths.push_back(read.Seq.length());
    }

    const float median = Median(&lengths);
    size_t maxLen = 2 * static_cast<size_t>(median);

    // if it's too short, return nothing
    if (median < static_cast<float>(minLength))
        return results;

    for (const auto& read : reads)
    {
        if (read.Seq.length() < maxLen)
            results.push_back(&read);
    }

    return results;
}

template<typename TRead>
boost::optional<ConsensusCore::MappedArrowRead>
ExtractMappedRead(const TRead& read,
                  const PoaAlignmentSummary& summary,
                  const size_t minLength)
{
    using ConsensusCore::FORWARD_STRAND;
    using ConsensusCore::REVERSE_STRAND;

    const size_t tplStart  = summary.ExtentOnConsensus.Left();
    const size_t tplEnd    = summary.ExtentOnConsensus.Right();
    const size_t readStart = summary.ExtentOnRead.Left();
    const size_t readEnd   = summary.ExtentOnRead.Right();

    if (readStart > readEnd || readEnd - readStart < minLength)
    {
        PBLOG_DEBUG << "Skipping read " << read.Id
                    << ", too short (<" << minLength << ')';
        return boost::none;
    }

    ConsensusCore::ArrowSequenceFeatures features(
            read.Seq.substr(readStart, readEnd - readStart));

    ConsensusCore::MappedArrowRead mappedRead(
            ConsensusCore::ArrowRead(features, read.Id, "N/A"),
            summary.ReverseComplementedRead ? REVERSE_STRAND : FORWARD_STRAND,
            tplStart,
            tplEnd);

    return boost::make_optional(mappedRead);
}

inline
std::string QVsToASCII(const std::vector<int>& qvs)
{
    std::string res;

    for (const int qv : qvs)
    {
        res.push_back(static_cast<char>(std::min(std::max(0, qv), 93) + 33));
    }

    return res;
}

#if 0
template<typename TRead>
bool ReadAccuracyDescending(const std::pair<size_t, const TRead*>& a,
                            const std::pair<size_t, const TRead*>& b)
{
    return a.second->ReadAccuracy > b.second->ReadAccuracy;
}
#endif

} // namespace anonymous


template<typename TRead>
std::string PoaConsensus(const std::vector<const TRead*>& reads,
                         std::vector<SparsePoa::ReadKey>* readKeys,
                         std::vector<PoaAlignmentSummary>* summaries,
                         const size_t maxPoaCov)
{
    SparsePoa poa;
    size_t cov = 0;

#if 0
    // create a vector of indices into the original reads vector,
    //   sorted by the ReadAccuracy in descending order
    std::vector<std::pair<size_t, const TRead*>> sorted;

    for (size_t i = 0; i < reads.size(); ++i)
        sorted.emplace_back(std::make_pair(i, reads[i]));

    std::sort(sorted.begin(), sorted.end(), ReadAccuracyDescending<TRead>);
#endif

    // initialize readKeys and resize
    readKeys->clear();
    // readKeys->resize(sorted.size());

    for (const auto read : reads)
    {
        SparsePoa::ReadKey key = poa.OrientAndAddRead(read->Seq);
        // SparsePoa::ReadKey key = poa.OrientAndAddRead(read.second->Seq);
        // readKeys->at(read.first) = key;
        readKeys->emplace_back(key);
        if (key >= 0 && (++cov) >= maxPoaCov)
            break;
    }

    // at least 50% of the reads should cover
    // TODO(lhepler) revisit this minimum coverage equation
    const size_t minCov = (cov < 5) ? 1 : (cov + 1) / 2 - 1;
    return poa.FindConsensus(minCov, &(*summaries))->Sequence;
}


// pass unique_ptr by reference to satisfy finickyness wrt move semantics in <future>
//   but then take ownership here with a local unique_ptr
template<typename TChunk, typename TResult>
ResultType<TResult> Consensus(std::unique_ptr<std::vector<TChunk>>& chunksRef,
                              const ConsensusSettings& settings)
{
    using namespace ConsensusCore::Arrow;

    auto chunks(std::move(chunksRef));
    ResultType<TResult> results;

    if (!chunks)
        return results;

    for (const auto& chunk : *chunks)
    {
        try
        {
            Timer timer;
            auto reads = FilterReads(chunk.Reads, settings.MinLength);

            if (reads.empty())
            {
                results.NoSubreads += 1;
                PBLOG_DEBUG << "Skipping " << chunk.Id
                            << ", no high quality subreads available";
                continue;
            }

            std::vector<SparsePoa::ReadKey> readKeys;
            std::vector<PoaAlignmentSummary> summaries;
            std::string poaConsensus = PoaConsensus(reads, &readKeys, &summaries,
                                                    settings.MaxPoaCoverage);

            if (poaConsensus.length() < settings.MinLength)
            {
                results.TooShort += 1;
                PBLOG_DEBUG << "Skipping " << chunk.Id
                            << ", initial consensus too short (<"
                            << settings.MinLength << ')';
                continue;
            }

            // setup the arrow scorer
            ContextParameters ctxParams(chunk.SignalToNoise);
            ArrowConfig config(ctxParams, BandingOptions(12.5));
            ArrowMultiReadMutationScorer scorer(config, poaConsensus);
            std::vector<int32_t> statusCounts(OTHER + 1, 0);
            const size_t nReads = readKeys.size();
            size_t nPasses = 0, nDropped = 0;

            // add the reads to the scorer
            for (size_t i = 0; i < nReads; ++i)
            {
                // skip unadded reads
                if (readKeys[i] < 0)
                    continue;

                if (auto mr = ExtractMappedRead(*reads[i], summaries[readKeys[i]], settings.MinLength))
                {
                    auto status = scorer.AddRead(*mr, settings.MinZScore);

                    // increment the status count
                    statusCounts[status] += 1;

                    if (status == SUCCESS &&
                        reads[i]->Flags & BAM::ADAPTER_BEFORE &&
                        reads[i]->Flags & BAM::ADAPTER_AFTER)
                    {
                        ++nPasses;
                    }
                    else if (status != SUCCESS)
                    {
                        ++nDropped;
                        PBLOG_DEBUG << "Skipping read " << mr->Name
                                    << ", " << AddReadResultNames[status];
                    }
                }
            }

            if (nPasses < settings.MinPasses)
            {
                results.TooFewPasses += 1;
                PBLOG_DEBUG << "Skipping " << chunk.Id
                            << ", insufficient number of passes ("
                            << nPasses << '<' << settings.MinPasses << ')';
                continue;
            }

            const double fracDropped = static_cast<double>(nDropped) / nReads;
            if (fracDropped > settings.MaxDropFraction)
            {
                results.TooManyUnusable += 1;
                PBLOG_DEBUG << "Skipping " << chunk.Id
                            << ", too high a fraction of unusable subreads ("
                            << fracDropped << '>' << settings.MaxDropFraction << ')';
                continue;
            }

            // get the original zscores
            const auto zdata = scorer.ZScores();

            // find consensus!!
            size_t nTested = 0, nApplied = 0;
            if (!RefineConsensus(scorer, &nTested, &nApplied))
            {
                results.NonConvergent += 1;
                PBLOG_DEBUG << "Skipping " << chunk.Id
                            << ", failed to converge";
                continue;
            }

            // compute predicted accuracy
            double predAcc = 0.0;
            std::vector<int> qvs = ConsensusQVs(scorer);
            for (const int qv : qvs)
            {
                predAcc += pow(10.0, static_cast<double>(qv) / -10.0);
            }
            predAcc = 1.0 - predAcc / qvs.size();

            if (predAcc < settings.MinPredictedAccuracy)
            {
                results.PoorQuality += 1;
                PBLOG_DEBUG << "Skipping " << chunk.Id
                            << ", failed to meet minimum predicted accuracy ("
                            << predAcc << '<' << settings.MinPredictedAccuracy << ')';
                continue;
            }

            // return resulting sequence!!
            results.Success += 1;
            results.emplace_back(
                TResult
                {
                    chunk.Id,
                    scorer.Template(),
                    QVsToASCII(qvs),
                    nPasses,
                    predAcc,
                    zdata.first.first,
                    zdata.first.second,
                    zdata.second,
                    statusCounts,
                    nTested,
                    nApplied,
                    chunk.SignalToNoise,
                    timer.ElapsedMilliseconds()
                });
        }
        catch (...)
        {
            results.Other += 1;
            PBLOG_ERROR << "Skipping " << chunk.Id
                        << ", caught exception during processing";
        }
    }

    return results;
}

} // namespace CCS
} // namespace PacBio
