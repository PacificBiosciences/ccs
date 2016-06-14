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
#include <functional>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <boost/optional.hpp>

#include <OptionParser.h>

#include <pacbio/consensus/Integrator.h>
#include <pacbio/consensus/Polish.h>
#include <pacbio/consensus/poa/PoaConsensus.h>

#include <pbbam/Accuracy.h>
#include <pbbam/LocalContextFlags.h>

#include <pacbio/ccs/Logging.h>
#include <pacbio/ccs/ReadId.h>
#include <pacbio/ccs/SparsePoa.h>
#include <pacbio/ccs/SubreadResultCounter.h>
#include <pacbio/ccs/Timer.h>

namespace PacBio {
namespace CCS {

using SNR = PacBio::Consensus::SNR;
using QualityValues = PacBio::Consensus::QualityValues;
using LocalContextFlags = PacBio::BAM::LocalContextFlags;
using Accuracy = PacBio::BAM::Accuracy;
using StrandEnum = PacBio::Consensus::StrandEnum;

namespace OptionNames {
// constexpr auto MaxPoaCoverage       = "maxPoaCoverage";
constexpr auto MaxLength = "maxLength";
constexpr auto MinLength = "minLength";
constexpr auto MinPasses = "minPasses";
constexpr auto MinPredictedAccuracy = "minPredictedAccuracy";
constexpr auto MinZScore = "minZScore";
constexpr auto MaxDropFraction = "maxDropFraction";
constexpr auto NoPolish = "noPolish";
constexpr auto MinReadScore = "minReadScore";
constexpr auto MinSnr = "minSnr";
constexpr auto ByStrand = "byStrand";
}  // namespace OptionNames

struct ConsensusSettings
{
    size_t MaxPoaCoverage;
    size_t MaxLength;
    size_t MinLength;
    size_t MinPasses;
    double MinPredictedAccuracy;
    double MinZScore;
    double MaxDropFraction;
    bool ByStrand;
    bool NoPolish;
    double MinReadScore;
    double MinSNR;

    ConsensusSettings(const optparse::Values& options);

    static void AddOptions(optparse::OptionParser* const parser)
    {
        const std::string em = "--";
        // TODO(lhepler) implement alignment to POA
        // parser->add_option(em +
        // OptionNames::MaxPoaCoverage).type("int").set_default(1024).help("Maximum number of
        // subreads to use when building POA. Default = %default");
        parser->add_option(em + OptionNames::MaxLength)
            .type("int")
            .set_default(7000)
            .help("Maximum length of subreads to use for generating CCS. Default = %default");
        parser->add_option(em + OptionNames::MinLength)
            .type("int")
            .set_default(10)
            .help("Minimum length of subreads to use for generating CCS. Default = %default");
        parser->add_option(em + OptionNames::MinPasses)
            .type("int")
            .set_default(3)
            .help("Minimum number of subreads required to generate CCS. Default = %default");
        parser->add_option(em + OptionNames::MinPredictedAccuracy)
            .type("float")
            .set_default(0.90)
            .help("Minimum predicted accuracy in [0, 1]. Default = %default");
        parser->add_option(em + OptionNames::MinZScore)
            .type("float")
            .set_default(-3.5)
            .help("Minimum z-score to use a subread. NaN disables this filter. Default = %default");
        parser->add_option(em + OptionNames::MaxDropFraction)
            .type("float")
            .set_default(0.34)
            .help(
                "Maximum fraction of subreads that can be dropped before giving up. Default = "
                "%default");
        parser->add_option(em + OptionNames::MinSnr)
            .type("float")
            .set_default(
                3.75)  // See https://github.com/PacificBiosciences/pbccs/issues/86 for a more
                       // detailed discussion of this default.
            .help("Minimum SNR of input subreads. Default = %default");
        parser->add_option(em + OptionNames::MinReadScore)
            .type("float")
            .set_default(0.75)
            .help("Minimum read score of input subreads. Default = %default");
        parser->add_option(em + OptionNames::ByStrand)
            .action("store_true")
            .help("Generate a consensus for each strand. Default = false");
        parser->add_option(em + OptionNames::NoPolish)
            .action("store_true")
            .help(
                "Only output the initial template derived from the POA (faster, less accurate). "
                "Default = false");
    }
};

template <typename TId>
struct ReadType
{
    TId Id;
    std::string Seq;
    std::vector<uint8_t> IPD;
    std::vector<uint8_t> PulseWidth;
    LocalContextFlags Flags;
    Accuracy ReadAccuracy;
    // TODO (move SNR and Chemistry here, eventually)
    // SNR SignalToNoise;
    // std::string Chemistry;
};

template <typename TId, typename TRead>
struct ChunkType
{
    TId Id;
    std::vector<TRead> Reads;
    SNR SignalToNoise;
    std::string Chemistry;
    boost::optional<std::tuple<int16_t, int16_t, uint8_t>> Barcodes;
};

struct ConsensusType
{
    ReadId Id;
    boost::optional<StrandEnum> Strand;
    std::string Sequence;
    QualityValues QVs;
    size_t NumPasses;
    double PredictedAccuracy;
    double AvgZScore;
    std::vector<double> ZScores;
    std::vector<int32_t> StatusCounts;
    size_t MutationsTested;
    size_t MutationsApplied;
    SNR SignalToNoise;
    float ElapsedMilliseconds;
    boost::optional<std::tuple<int16_t, int16_t, uint8_t>> Barcodes;
};

template <typename TConsensus>
class ResultType : public std::vector<TConsensus>
{
public:
    size_t Success;
    size_t PoorSNR;
    size_t NoSubreads;
    size_t TooLong;
    size_t TooShort;
    size_t TooFewPasses;
    size_t TooManyUnusable;
    size_t NonConvergent;
    size_t PoorQuality;
    size_t ExceptionThrown;
    SubreadResultCounter SubreadCounter;

    ResultType()
        : Success{0}
        , PoorSNR{0}
        , NoSubreads{0}
        , TooLong{0}
        , TooShort{0}
        , TooFewPasses{0}
        , TooManyUnusable{0}
        , NonConvergent{0}
        , PoorQuality{0}
        , ExceptionThrown{0}
        , SubreadCounter{}
    {
    }

    ResultType<TConsensus>& operator+=(const ResultType<TConsensus>& other)
    {
        Success += other.Success;
        PoorSNR += other.PoorSNR;
        NoSubreads += other.NoSubreads;
        TooLong += other.TooLong;
        TooShort += other.TooShort;
        TooManyUnusable += other.TooManyUnusable;
        TooFewPasses += other.TooFewPasses;
        NonConvergent += other.NonConvergent;
        PoorQuality += other.PoorQuality;
        ExceptionThrown += other.ExceptionThrown;
        SubreadCounter += other.SubreadCounter;
        return *this;
    }

    size_t Total() const
    {
        return (Success + PoorSNR + NoSubreads + TooShort + TooManyUnusable + TooFewPasses +
                NonConvergent + PoorQuality + ExceptionThrown);
    }
};

namespace {  // anonymous

template <typename T>
float Median(std::vector<T>* vs)
{
    const size_t n = vs->size();
    std::sort(vs->begin(), vs->end());
    if (n % 2 == 1) return static_cast<float>(vs->at(n / 2));
    return 0.5 * (vs->at(n / 2 - 1) + vs->at(n / 2));
}

template <typename TRead>
std::vector<const TRead*> FilterReads(const std::vector<TRead>& reads,
                                      const ConsensusSettings& settings,
                                      SubreadResultCounter* resultCounter)
{
    // This is a count of subreads removed for bing too short, or too long.
    std::vector<const TRead*> results;

    if (reads.empty()) return results;

    std::vector<size_t> lengths;
    size_t longest = 0;

    // get the lengths for all full-length subreads
    for (const auto& read : reads) {
        longest = std::max(longest, read.Seq.length());
        if ((read.Flags & BAM::ADAPTER_BEFORE) && (read.Flags & BAM::ADAPTER_AFTER) &&
            read.ReadAccuracy >= settings.MinReadScore)
            lengths.emplace_back(read.Seq.length());
    }

    // nonexistent median is just the greatest observed length
    const float median = lengths.empty() ? static_cast<float>(longest) : Median(&lengths);
    size_t maxLen = std::min(2 * static_cast<size_t>(median), settings.MaxLength);

    // if it's too short, return nothing
    if (median < static_cast<float>(settings.MinLength)) {
        resultCounter->FilteredBySize += reads.size();
        return results;
    }
    results.reserve(reads.size());

    for (const auto& read : reads) {
        // if the median exists, then this filters stuff,
        //   otherwise it's twice the longest read and is always true

        if (read.ReadAccuracy < settings.MinReadScore) {
            resultCounter->BelowMinQual += 1;
            results.emplace_back(nullptr);
        } else if (read.Seq.length() < maxLen) {
            results.emplace_back(&read);
        } else {
            resultCounter->FilteredBySize += 1;
            results.emplace_back(nullptr);
        }
    }

    // TODO(lhepler): incorporate per-subread quality here
    // End-to-end reads take priority, hence the lexicographical sort;
    //   always take the read with the least deviation from the median.
    //   In the case of no median, longer reads are prioritized.
    const auto lexForm = [median](const TRead* read) {
        const float l = static_cast<float>(read->Seq.length());
        const float v = std::min(l / median, median / l);

        if (read->Flags & BAM::ADAPTER_BEFORE && read->Flags & BAM::ADAPTER_AFTER)
            return std::make_tuple(v, 0.0f);

        return std::make_tuple(0.0f, v);
    };

    std::stable_sort(results.begin(), results.end(),
                     [&lexForm](const TRead* lhs, const TRead* rhs) {
                         if (lhs == nullptr)
                             return false;
                         else if (rhs == nullptr)
                             return true;

                         return lexForm(lhs) > lexForm(rhs);
                     });

    return results;
}

template <typename TRead>
boost::optional<PacBio::Consensus::MappedRead> ExtractMappedRead(
    const TRead& read, const PacBio::Consensus::SNR& snr, const std::string& chem,
    const PoaAlignmentSummary& summary, const size_t poaLength, const ConsensusSettings& settings,
    SubreadResultCounter* resultCounter)
{
    constexpr size_t kStickyEnds = 7;

    size_t readStart = summary.ExtentOnRead.Left();
    size_t readEnd = summary.ExtentOnRead.Right();
    size_t tplStart = summary.ExtentOnConsensus.Left();
    size_t tplEnd = summary.ExtentOnConsensus.Right();

    // if we're ADAPTER_BEFORE and _AFTER and mapped nearly end-to-end,
    //   just make it end to end (but for each side, respectively)
    if (summary.ReverseComplementedRead) {
        if (read.Flags & BAM::ADAPTER_BEFORE && (poaLength - tplEnd) <= kStickyEnds)
            tplEnd = poaLength;
        if (read.Flags & BAM::ADAPTER_AFTER && tplStart <= kStickyEnds) tplStart = 0;
    } else {
        if (read.Flags & BAM::ADAPTER_BEFORE && tplStart <= kStickyEnds) tplStart = 0;
        if (read.Flags & BAM::ADAPTER_AFTER && (poaLength - tplEnd) <= kStickyEnds)
            tplEnd = poaLength;
    }

    if (readStart > readEnd || readEnd - readStart < settings.MinLength) {
        resultCounter->FilteredBySize += 1;
        PBLOG_DEBUG << "Skipping read " << read.Id << ", too short (<" << settings.MinLength << ')';
        return boost::none;
    } else if (readEnd - readStart > settings.MaxLength) {
        resultCounter->FilteredBySize += 1;
        PBLOG_DEBUG << "Skipping read " << read.Id << ", too long (>" << settings.MaxLength << ')';
        return boost::none;
    }

    PacBio::Consensus::MappedRead mappedRead(
        PacBio::Consensus::Read(
            read.Id, read.Seq.substr(readStart, readEnd - readStart),
            std::vector<uint8_t>(read.IPD.begin() + readStart, read.IPD.begin() + readEnd),
            std::vector<uint8_t>(read.PulseWidth.begin() + readStart,
                                 read.PulseWidth.begin() + readEnd),
            snr, chem),
        summary.ReverseComplementedRead ? StrandEnum::REVERSE : StrandEnum::FORWARD, tplStart,
        tplEnd, (tplStart == 0) ? true : false, (tplEnd == poaLength) ? true : false);

    return boost::make_optional(mappedRead);
}

#if 0
template<typename TRead>
bool ReadAccuracyDescending(const std::pair<size_t, const TRead*>& a,
                            const std::pair<size_t, const TRead*>& b)
{
    return a.second->ReadAccuracy > b.second->ReadAccuracy;
}
#endif

}  // namespace anonymous

template <typename TRead>
std::string PoaConsensus(const std::vector<const TRead*>& reads,
                         std::vector<SparsePoa::ReadKey>* readKeys,
                         std::vector<PoaAlignmentSummary>* summaries, const size_t maxPoaCov)
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

    for (const auto read : reads) {
        SparsePoa::ReadKey key = (read == nullptr) ? -1 : poa.OrientAndAddRead(read->Seq);
        // SparsePoa::ReadKey key = poa.OrientAndAddRead(read.second->Seq);
        // readKeys->at(read.first) = key;
        readKeys->emplace_back(key);
        if (key >= 0 && (++cov) >= maxPoaCov) break;
    }

    // at least 50% of the reads should cover
    // TODO(lhepler) revisit this minimum coverage equation
    const size_t minCov = (cov < 5) ? 1 : (cov + 1) / 2 - 1;
    return poa.FindConsensus(minCov, &(*summaries))->Sequence;
}

// pass unique_ptr by reference to satisfy finickyness wrt move semantics in <future>
//   but then take ownership here with a local unique_ptr
template <typename TChunk>
ResultType<ConsensusType> Consensus(std::unique_ptr<std::vector<TChunk>>& chunksRef,
                                    const ConsensusSettings& settings)
{
    using namespace PacBio::Consensus;

    auto chunks(std::move(chunksRef));
    ResultType<ConsensusType> result;

    if (!chunks) return result;
    // We should only be dealing with chunks of size 1
    if (chunks->size() != 1) {
        throw std::runtime_error("CCS chunk was of size != 1");
    }
    const auto& chunk = chunks->at(0);

    try {
        Timer timer;

        // Do read level SNR filtering first
        auto& snr = chunk.SignalToNoise;
        auto minSNR = std::min(std::min(snr.A, snr.C), std::min(snr.G, snr.T));
        if (minSNR < settings.MinSNR) {
            result.SubreadCounter.ZMWBelowMinSNR += chunk.Reads.size();
            result.PoorSNR += 1;
            return result;
        }

        auto reads = FilterReads(chunk.Reads, settings, &result.SubreadCounter);

        if (reads.empty() ||  // Check if subread are present
            std::accumulate(reads.begin(), reads.end(), 0, std::plus<bool>()) == 0) {
            result.NoSubreads += 1;
            PBLOG_DEBUG << "Skipping " << chunk.Id << ", no high quality subreads available";
            return result;
        }

        // If it is not possible to exceed the minPasses requirement, we will bail here before
        //   generating the POA, filling the matrices and performing all the other checks
        size_t possiblePasses = 0;
        size_t activeReads = 0;
        for (size_t i = 0; i < reads.size(); ++i) {
            if (reads[i] != nullptr) {
                activeReads += 1;
                if (reads[i]->Flags & BAM::ADAPTER_BEFORE && reads[i]->Flags & BAM::ADAPTER_AFTER) {
                    possiblePasses += 1;
                }
            }
        }

        if (possiblePasses < settings.MinPasses) {
            result.TooFewPasses += 1;
            result.SubreadCounter.ZMWNotEnoughSubReads += activeReads;
            PBLOG_DEBUG << "Skipping " << chunk.Id << ", not enough possible passes ("
                        << possiblePasses << '<' << settings.MinPasses << ')';
            return result;
        }

        std::vector<SparsePoa::ReadKey> readKeys;
        std::vector<PoaAlignmentSummary> summaries;
        std::string poaConsensus =
            PoaConsensus(reads, &readKeys, &summaries, settings.MaxPoaCoverage);

        if (poaConsensus.length() < settings.MinLength) {
            result.TooShort += 1;
            result.SubreadCounter.Other += activeReads;
            PBLOG_DEBUG << "Skipping " << chunk.Id << ", initial consensus too short (<"
                        << settings.MinLength << ')';
        } else if (poaConsensus.length() > settings.MaxLength) {
            result.TooLong += 1;
            result.SubreadCounter.Other += activeReads;
            PBLOG_DEBUG << "Skipping " << chunk.Id << ", initial consensus too long (>"
                        << settings.MaxLength << ')';
        } else {
            if (settings.NoPolish) {
                const size_t len = poaConsensus.length();
                // generate dummy QVs, will use 20
                // TODO(lhepler): should we use different values for delQVs, insQVs, and subQVs?
                QualityValues qvs{std::vector<int>(len, 20), std::vector<int>(len, 20),
                                  std::vector<int>(len, 20), std::vector<int>(len, 20)};
                result.Success += 1;
                result.SubreadCounter.Success += activeReads;
                result.emplace_back(ConsensusType{
                    chunk.Id, boost::none, poaConsensus, qvs, possiblePasses, 0, 0,
                    std::vector<double>(1), result.SubreadCounter.ReturnCountsAsArray(), 0, 0,
                    chunk.SignalToNoise, timer.ElapsedMilliseconds(), chunk.Barcodes});
            } else {
                const auto mkConsensus = [&](const boost::optional<StrandEnum> strand) {
                    // give this consensus attempt a name we can refer to
                    std::string chunkName(chunk.Id);
                    if (strand && *strand == StrandEnum::FORWARD) chunkName += " [fwd]";
                    if (strand && *strand == StrandEnum::REVERSE) chunkName += " [rev]";

                    try {
                        // setup the arrow integrator
                        IntegratorConfig cfg(settings.MinZScore, 12.5);
                        MonoMolecularIntegrator ai(poaConsensus, cfg, chunk.SignalToNoise,
                                                   chunk.Chemistry);
                        const size_t nReads = readKeys.size();
                        size_t nPasses = 0, nDropped = 0;

                        // If this ZMW could possibly pass,  add the reads to the integrator
                        for (size_t i = 0; i < nReads; ++i) {
                            // skip unadded reads
                            if (readKeys[i] < 0) continue;

                            if (auto mr = ExtractMappedRead(*reads[i], chunk.SignalToNoise,
                                                            chunk.Chemistry, summaries[readKeys[i]],
                                                            poaConsensus.length(), settings,
                                                            &result.SubreadCounter)) {
                                // skip reads not belonging to this strand, if we're --byStrand
                                if (strand && mr->Strand != *strand) continue;
                                auto status = ai.AddRead(*mr);
                                // increment the status count
                                result.SubreadCounter.AddResult(status);
                                if (status == AddReadResult::SUCCESS &&
                                    reads[i]->Flags & BAM::ADAPTER_BEFORE &&
                                    reads[i]->Flags & BAM::ADAPTER_AFTER) {
                                    nPasses += 1;
                                } else if (status != AddReadResult::SUCCESS) {
                                    nDropped += 1;
                                    PBLOG_DEBUG << "Skipping read " << mr->Name << ", " << status;
                                }
                            }
                        }

                        if (nPasses < settings.MinPasses) {
                            // Reassign all the successful reads to the other category
                            result.SubreadCounter.AssignSuccessToOther();
                            result.TooFewPasses += 1;
                            PBLOG_DEBUG << "Skipping " << chunkName
                                        << ", insufficient number of passes (" << nPasses << '<'
                                        << settings.MinPasses << ')';
                            return;
                        }

                        // this is hairy, but also relatively straightforward, so bear with me:
                        //   if we're not doing strand-specific consensus, the total number of
                        //   available reads is just nReads. If we're doing byStrand though,
                        //   then the number of available reads is those that mapped to this
                        //   strand, plus half of those that didn't map to the POA (we assume).
                        const size_t nAvail =
                            (!strand) ? nReads
                                      : (std::count_if(
                                             readKeys.cbegin(), readKeys.cend(),
                                             [&](const SparsePoa::ReadKey key) {
                                                 return key >= 0 &&
                                                        summaries[key].ReverseComplementedRead ==
                                                            (*strand == StrandEnum::REVERSE);
                                             }) +
                                         std::count_if(
                                             readKeys.cbegin(), readKeys.cend(),
                                             [](const SparsePoa::ReadKey key) { return key < 0; }) /
                                             2);

                        const double fracDropped = static_cast<double>(nDropped) / nAvail;
                        if (fracDropped > settings.MaxDropFraction) {
                            result.TooManyUnusable += 1;
                            result.SubreadCounter.AssignSuccessToOther();
                            PBLOG_DEBUG << "Skipping " << chunkName
                                        << ", too high a fraction of unusable subreads ("
                                        << fracDropped << '>' << settings.MaxDropFraction << ')';
                            return;
                        }

                        const double zAvg = ai.AvgZScore();
                        const auto zScores = ai.ZScores();

                        // find consensus!!
                        size_t nTested, nApplied;
                        bool polished;
                        std::tie(polished, nTested, nApplied) = Polish(&ai, PolishConfig());

                        if (!polished) {
                            result.NonConvergent += 1;
                            result.SubreadCounter.AssignSuccessToOther();
                            PBLOG_DEBUG << "Skipping " << chunkName << ", failed to converge";
                            return;
                        }

                        // compute predicted accuracy
                        double predAcc = 0.0;
                        QualityValues qvs = ConsensusQVs(ai);
                        for (const int qv : qvs.Qualities) {
                            predAcc += pow(10.0, static_cast<double>(qv) / -10.0);
                        }
                        predAcc = 1.0 - predAcc / qvs.Qualities.size();

                        if (predAcc < settings.MinPredictedAccuracy) {
                            result.PoorQuality += 1;
                            result.SubreadCounter.AssignSuccessToOther();
                            PBLOG_DEBUG << "Skipping " << chunkName
                                        << ", failed to meet minimum predicted accuracy ("
                                        << predAcc << '<' << settings.MinPredictedAccuracy << ')';
                            return;
                        }

                        // return resulting sequence!!
                        result.Success += 1;
                        result.emplace_back(ConsensusType{
                            chunk.Id, strand, std::string(ai), std::move(qvs), nPasses, predAcc,
                            zAvg, zScores, result.SubreadCounter.ReturnCountsAsArray(), nTested,
                            nApplied, chunk.SignalToNoise, timer.ElapsedMilliseconds(),
                            chunk.Barcodes});
                    } catch (const std::exception& e) {
                        result.ExceptionThrown += 1;
                        PBLOG_ERROR << "Skipping " << chunkName << ", caught exception: '"
                                    << e.what() << "\'";
                    }
                };

                if (settings.ByStrand) {
                    mkConsensus(StrandEnum::FORWARD);
                    mkConsensus(StrandEnum::REVERSE);
                } else {
                    mkConsensus(boost::none);
                }
            }
        }
    } catch (const std::exception& e) {
        result.ExceptionThrown += 1;
        PBLOG_ERROR << "Skipping " << chunk.Id << ", caught exception: '" << e.what() << "\'";
    } catch (...) {
        // This should NEVER happen. Only here as a guard, if this is ever printed someone
        // goofed
        // up by throwing something that didn't derive from std::exception.
        result.ExceptionThrown += 1;
        PBLOG_ERROR << "Skipping " << chunk.Id << ", caught unknown exception type";
    }

    return result;
}

}  // namespace CCS
}  // namespace PacBio
