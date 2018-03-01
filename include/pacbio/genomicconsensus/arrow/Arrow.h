// Author: Derek Barnett

#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <pbbam/BamRecord.h>
#include <boost/optional.hpp>

#include <pacbio/align/AffineAlignment.h>
#include <pacbio/align/LinearAlignment.h>
#include <pacbio/align/PairwiseAlignment.h>
#include <pacbio/consensus/Integrator.h>
#include <pacbio/consensus/Polish.h>
#include <pacbio/data/Interval.h>
#include <pacbio/data/Read.h>
#include <pacbio/data/Sequence.h>
#include <pacbio/data/State.h>
#include <pacbio/data/internal/BaseEncoding.h>
#include <pacbio/denovo/PoaConsensus.h>
#include <pacbio/genomicconsensus/Consensus.h>
#include <pacbio/genomicconsensus/Input.h>
#include <pacbio/genomicconsensus/ReferenceWindow.h>
#include <pacbio/genomicconsensus/Settings.h>
#include <pacbio/genomicconsensus/Variant.h>
#include <pacbio/genomicconsensus/WindowResult.h>

namespace PacBio {
namespace GenomicConsensus {

struct Arrow
{
    static size_t Clamp(const size_t pos, const size_t min, const size_t max)
    {
        if (pos < min) return min;
        if (pos > max) return max;
        return pos;
    };

    static size_t Median(std::vector<size_t> v)
    {
        const auto mid = v.size() / 2;
        auto nth = v.begin() + mid;
        std::nth_element(v.begin(), nth, v.end());
        if (v.size() % 2 == 0) {
            // avg of middle two
            const auto first = *nth;
            std::nth_element(v.begin(), --nth, v.end());
            const auto second = *nth;
            return (first + second) / 2;
        }
        return *nth;
    }

    static std::unique_ptr<const PacBio::Poa::PoaConsensus> MakePoaConsensus(
        std::vector<std::string>&& fwdSequences, const Settings& settings)
    {
        using AlignMode = PacBio::Align::AlignMode;
        using PoaConsensus = PacBio::Poa::PoaConsensus;

        std::vector<size_t> seqLengths;
        seqLengths.reserve(fwdSequences.size());
        for (const auto& seq : fwdSequences)
            seqLengths.push_back(seq.size());

        const auto median = Median(std::move(seqLengths));
        std::vector<std::string> ordSeqs = std::move(fwdSequences);
        std::sort(ordSeqs.begin(), ordSeqs.end(),
                  [median](const std::string& lhs, const std::string& rhs) {
                      const auto lhsSize = static_cast<int>(lhs.size());
                      const auto rhsSize = static_cast<int>(rhs.size());
                      const auto intMedian = static_cast<int>(median);
                      const auto diffLeft = std::abs(lhsSize - intMedian);
                      const auto diffRight = std::abs(rhsSize - intMedian);
                      return diffLeft < diffRight;
                  });
        ordSeqs.resize(std::distance(ordSeqs.begin(), ordSeqs.begin() + settings.maxPoaCoverage));

        const auto poaConfig = PacBio::Poa::DefaultPoaConfig(AlignMode::GLOBAL);
        const auto cov = ordSeqs.size();
        const auto minCov = (cov < 5 ? 1 : ((cov + 1) / 2 - 1));
        return std::unique_ptr<const PacBio::Poa::PoaConsensus>{
            PoaConsensus::FindConsensus(ordSeqs, poaConfig, minCov)};
    }

    static std::vector<std::string> FilteredForwardSequences(
        const std::vector<PacBio::BAM::BamRecord>& reads, const ReferenceWindow& window)
    {
        using BamRecord = PacBio::BAM::BamRecord;
        using Orientation = PacBio::BAM::Orientation;
        using Strand = PacBio::BAM::Strand;

        const auto spansReferenceRange = [](const BamRecord& read, const ReferenceWindow& window) {
            const auto tStart = static_cast<size_t>(read.ReferenceStart());
            const auto tEnd = static_cast<size_t>(read.ReferenceEnd());
            assert(window.Start() <= window.End());
            return (tStart <= window.Start() && tEnd >= window.End());
        };

        std::vector<std::string> result;
        for (const auto& read : reads) {
            if (read.AlignedStrand() == Strand::FORWARD && spansReferenceRange(read, window))
                result.push_back(read.Sequence(Orientation::NATIVE, false));  // excise soft clips?
        }
        return result;
    }

    static std::vector<PacBio::Data::MappedRead> MakeMappedReads(
        const std::vector<PacBio::BAM::BamRecord>& reads, const ReferenceWindow& window)
    {
        using MappedRead = PacBio::Data::MappedRead;
        using SNR = PacBio::Data::SNR;
        using Read = PacBio::Data::Read;
        using Strand = PacBio::BAM::Strand;
        using StrandType = PacBio::Data::StrandType;

        std::vector<MappedRead> mappedReads;
        mappedReads.reserve(reads.size());
        for (const auto& r : reads) {
            const auto seq = r.Sequence();
            std::vector<uint8_t> ipd;
            if (r.HasIPD())
                ipd = r.IPD().Encode();
            else
                ipd = std::vector<uint8_t>(seq.length(), 0);

            std::vector<uint8_t> pw;
            if (r.HasPulseWidth())
                pw = r.PulseWidth().Encode();
            else
                pw = std::vector<uint8_t>(seq.length(), 0);

            const auto bamSNR = r.SignalToNoise();
            const SNR snr(bamSNR[0], bamSNR[1], bamSNR[2], bamSNR[3]);
            mappedReads.emplace_back(
                Read{r.FullName(), seq, ipd, pw, snr, r.ReadGroup().SequencingChemistry()},
                (r.AlignedStrand() == Strand::FORWARD ? StrandType::FORWARD : StrandType::REVERSE),
                r.ReferenceStart() - window.Start(), r.ReferenceEnd() - window.Start());
        }
        return mappedReads;
    }

    static void LiftReads(const std::vector<int>& queryPositions,
                          std::vector<PacBio::Data::MappedRead>* mappedReads)
    {
        for (auto& read : *mappedReads) {
            read.TemplateStart = queryPositions[read.TemplateStart];
            read.TemplateEnd = queryPositions[read.TemplateEnd];
        }
    }

    static bool IsSufficientlyAccurate(const PacBio::Data::MappedRead& mr,
                                       const std::string& poaCss, const Settings& settings)
    {
        using PairwiseAlignment = PacBio::Align::PairwiseAlignment;
        using StrandType = PacBio::Data::StrandType;

        if (settings.minAccuracy <= 0.0) return true;

        const auto start = mr.TemplateStart;
        const auto end = mr.TemplateEnd;
        if (mr.Strand == StrandType::UNMAPPED) return false;

        auto tpl = poaCss.substr(start, end - start);
        if (mr.Strand == StrandType::REVERSE) tpl = PacBio::Data::ReverseComplement(tpl);

        std::unique_ptr<PairwiseAlignment> aln{PacBio::Align::AlignLinear(tpl, mr.Seq)};
        const auto transcript = aln->Transcript();

        // NOTE: this is *not* CIGAR 'M', but rather equivalent to CIGAR '='
        const size_t numErrors = std::count_if(transcript.begin(), transcript.end(),
                                               [](const char c) { return c != 'M'; });
        const size_t tLength = tpl.size();
        const double accuracy = 1.0 - (1.0 * std::min(numErrors, tLength) / tLength);
        return accuracy >= settings.minAccuracy;
    }

    static std::vector<uint8_t> ConsensusConfidence(PacBio::Consensus::Integrator& integrator)
    {
        //
        // Returns an array of QV values reflecting the consensus confidence
        // at each position specified.  If the `positions` argument is
        // omitted, confidence values are returned for all positions in the
        // consensus (str(ai)).
        //

        const auto cssQVs = PacBio::Consensus::ConsensusQualities(integrator);
        std::vector<uint8_t> confidence;
        confidence.reserve(cssQVs.size());
        for (const auto c : cssQVs)
            confidence.emplace_back(Clamp(c, 0, 93));
        return confidence;
    }

    static bool RefineConsensus(PacBio::Consensus::Integrator& integrator, const Settings& settings,
                                bool polishDiploid)
    {
        //
        // Given a MultiReadMutationScorer, identify and apply favorable
        // template mutations.  Return (consensus, didConverge) :: (str, bool)
        //

        using PolishConfig = PacBio::Consensus::PolishConfig;
        const auto config = PolishConfig{settings.maxIterations, settings.mutationSeparation,
                                         settings.mutationNeighborhood, polishDiploid};

        if (settings.maskRadius != 0) {
            PacBio::Consensus::Polish(&integrator, config);
            integrator.MaskIntervals(settings.maskRadius, settings.maskErrorRate);
        }
        const auto polishResult = PacBio::Consensus::Polish(&integrator, config);
        return polishResult.hasConverged;
    }

    static Consensus ConsensusForAlignments(const ReferenceWindow& window,
                                            const std::string& refSeq,
                                            const std::vector<PacBio::BAM::BamRecord>& reads,
                                            const Settings& settings,
                                            std::vector<PacBio::BAM::BamRecord>* readsUsed,
                                            const bool polishDiploid = true,  // check
                                            const bool polish = true)
    {
        //
        // Call consensus on this interval---without subdividing the interval
        // further.
        //
        // Returns a consensus object.
        //
        // Requires that clipping has already been done.
        //
        // If `draft` is provided, it will serve as the starting
        // point for polishing.  If not, the POA will be used to generate a
        // draft starting point.
        //
        // If `polish` is False, the arrow polishing procedure will not be
        // used, and the draft consensus will be returned.
        //
        // `alnsUsed` is an output parameter; if not None, it should be an
        // empty list on entry; on return from this function, the list will
        // contain the alns objects that were actually used to compute the
        // consensus (those not filtered out).
        //

        using Consensus = PacBio::GenomicConsensus::Consensus;
        using Integrator = PacBio::Consensus::Integrator;
        using IntegratorConfig = PacBio::Consensus::IntegratorConfig;
        using PairwiseAlignment = PacBio::Align::PairwiseAlignment;
        using State = PacBio::Data::State;
        using StrandType = PacBio::Data::StrandType;

        if (readsUsed) assert(readsUsed->empty());

        // Compute the POA consensus, which is our initial guess, and
        // should typically be > 99.5% accurate
        auto fwdSequences = FilteredForwardSequences(reads, window);
        assert(fwdSequences.size() >= settings.minPoaCoverage);

        std::string draft;
        try {
            const auto p = MakePoaConsensus(std::move(fwdSequences), settings);
            draft = p->Sequence;
        } catch (std::exception&) {
            //TODO: log
            return Consensus::NoCallConsensus(
                NoCallStyle::LOWERCASE_REFERENCE,  // double-check this
                window, refSeq);
        }

        // align draft to reference, map reads relative to POA consenus
        const std::unique_ptr<PairwiseAlignment> ga{PacBio::Align::Align(refSeq, draft)};
        auto mappedReads = MakeMappedReads(reads, window);
        const auto queryPositions = PacBio::Align::TargetToQueryPositions(*ga);
        LiftReads(queryPositions, &mappedReads);

        // Load the mapped reads into the mutation scorer, and iterate
        // until convergence.
        auto ai = Integrator{draft, IntegratorConfig{settings.minZScore}};
        size_t coverage = 0;
        for (size_t i = 0; i < mappedReads.size(); ++i) {
            const auto& mr = mappedReads.at(i);

            if ((mr.TemplateEnd <= mr.TemplateStart) || (mr.TemplateEnd - mr.TemplateStart < 2) ||
                (mr.Length() < 2)) {
                continue;
            }
            if (!IsSufficientlyAccurate(mr, draft, settings)) {
                auto tpl = draft.substr(mr.TemplateStart, mr.TemplateEnd - mr.TemplateStart);
                if (mr.Strand == StrandType::FORWARD)
                    ;
                else if (mr.Strand == StrandType::REVERSE)
                    tpl = PacBio::Data::ReverseComplement(tpl);
                else
                    tpl = "INACTIVE/UNMAPPED";
                //TODO: log
                continue;
            }

            if (ai.AddRead(mr) == State::VALID) {
                ++coverage;
                if (readsUsed) readsUsed->push_back(reads.at(i));
            }
        }

        if (coverage < settings.minPoaCoverage) {
            // log
            return Consensus::NoCallConsensus(NoCallStyle::LOWERCASE_REFERENCE, window, refSeq);
        }

        if (!polish) {
            return Consensus{
                window, draft, std::vector<uint8_t>(draft.size(), 0)
                // AI?
            };
        }

        // Iterate until covergence
        std::vector<uint8_t> confidence;
        std::string arrowCss;
        auto converged = RefineConsensus(ai, settings, false);
        if (converged) {
            arrowCss = static_cast<std::string>(ai);
            if (settings.computeConfidence)
                confidence = ConsensusConfidence(ai);
            else
                confidence = std::vector<uint8_t>(arrowCss.size(), 0);
        } else {
            //TODO: log
            return Consensus::NoCallConsensus(NoCallStyle::LOWERCASE_REFERENCE, window, refSeq);
        }

        if (settings.polishDiploid) {
            converged = RefineConsensus(ai, settings, true);
            if (converged) {
                arrowCss = static_cast<std::string>(ai);
                if (settings.computeConfidence)
                    confidence = ConsensusConfidence(ai);
                else
                    confidence = std::vector<uint8_t>(arrowCss.size(), 0);
            } else {
                //TODO: log
            }
        }

        return Consensus{
            window, arrowCss, confidence
            // AI?
        };
    }

    static std::vector<uint8_t> ProjectIntoRange(
        const std::vector<PacBio::Data::Interval>& intervals, const ReferenceWindow& window)
    {
        std::vector<uint8_t> result(window.Length(), 0);

        const auto winStart = window.Start();
        const auto winEnd = window.End();
        for (const auto& interval : intervals) {
            const size_t tStart = Clamp(interval.Left(), winStart, winEnd) - winStart;
            const size_t tEnd = Clamp(interval.Right(), winStart, winEnd) - winStart;
            for (size_t i = tStart; i < tEnd; ++i)
                result.at(i)++;
        }
        return result;
    }

    static std::vector<uint8_t> CoverageInWindow(const ReferenceWindow& window,
                                                 const std::vector<PacBio::BAM::BamRecord>& reads)
    {
        using Interval = PacBio::Data::Interval;

        std::vector<Interval> intervals;
        for (const auto& read : reads) {
            if (read.ReferenceName() == window.name) {
                const auto start = static_cast<size_t>(read.ReferenceStart());
                const auto end = static_cast<size_t>(read.ReferenceEnd());
                intervals.push_back({start, end});
            }
        }
        return ProjectIntoRange(intervals, window);
    }

    static std::vector<PacBio::Data::Interval> TranscriptIntervals(const std::string& transcript)
    {
        using Interval = PacBio::Data::Interval;

        std::vector<Interval> result;
        const auto tSize = transcript.size();
        if (tSize == 0) return result;
        assert(!transcript.empty());

        char previousChar = transcript.at(0);
        size_t currentRunStart = 0;
        size_t currentRunLength = 1;

        for (size_t i = 1; i < tSize; ++i) {
            const auto currentChar = transcript.at(i);
            if (currentChar == previousChar)
                ++currentRunLength;
            else {
                result.emplace_back(currentRunStart, currentRunStart + currentRunLength);
                currentRunStart = i;
                currentRunLength = 1;
            }
            previousChar = currentChar;
        }
        // push last interval & return
        result.emplace_back(currentRunStart, currentRunStart + currentRunLength);
        return result;
    }

    static UNANIMITY_CONSTEXPR const char* LookupIUPAC(const char c)
    {
        constexpr const std::array<const char*, 16> table{{"-", "A", "C", "AC", "G", "AG", "CG",
                                                           "ACG", "T", "AT", "CT", "ACT", "GT",
                                                           "AGT", "CGT", "ACGT"}};
        const auto ncbi4na = PacBio::Data::detail::NCBI4na::FromASCII(c);
        return table.at(ncbi4na.Data());
    }

    struct SplitupIUPACResult
    {
        std::string readSeq1;
        boost::optional<std::string> readSeq2 = boost::none;
        boost::optional<double> freq1 = boost::none;
        boost::optional<double> freq2 = boost::none;
    };

    static SplitupIUPACResult SplitupIUPAC(const std::string& css)
    {
        std::string listSeq1;
        std::string listSeq2;
        listSeq1.reserve(css.size());
        listSeq2.reserve(css.size());
        for (const auto c : css) {
            const std::string value{LookupIUPAC(c)};
            listSeq1.push_back(value.front());
            listSeq2.push_back(value.back());
        }

        SplitupIUPACResult result;
        if (listSeq1 == listSeq2) {
            // haploid
            result.readSeq1 = std::move(listSeq1);
        } else {
            // diploid
            result.readSeq1 = std::move(listSeq1);
            result.readSeq2 = std::move(listSeq2);
            result.freq1 = 0.5;
            result.freq2 = 0.5;
        }
        return result;
    };

    template <typename T>
    static bool FoundCode(const T& validCodes, const char code)
    {
        for (const char c : validCodes) {
            if (c == code) return true;
        }
        return false;
    }

    static std::vector<Variant> VariantsFromAlignment(
        const std::unique_ptr<PacBio::Align::PairwiseAlignment>& alignment,
        const ReferenceWindow& window, const boost::optional<std::vector<uint8_t>>& cssQvInWindow,
        const std::vector<uint8_t>& siteCoverage,
        const boost::optional<std::vector<uint8_t>>& effectiveSiteCoverage)
    {
        std::vector<Variant> variants;

        const auto refId = window.name;
        const auto refStart = window.Start();
        size_t refPos = refStart;
        size_t cssPos = 0;
        char refPrev = 'N';
        char cssPrev = 'N';

        constexpr const std::array<char, 5> validCodes{{'R', 'I', 'D', 'M', 'N'}};

        // We don't call variants where either the reference or css is 'N'
        const auto& target = alignment->Target();
        const auto& query = alignment->Query();
        auto transcript = alignment->Transcript();
        for (size_t i = 0; i < transcript.size(); ++i) {
            const char t = target.at(i);
            const char q = query.at(i);
            if (t == 'N' || q == 'N') transcript[i] = 'N';
        }

        boost::optional<Variant> variant;
        for (const auto& interval : TranscriptIntervals(transcript)) {

            const auto pos = interval.Left();
            const auto code = transcript.at(pos);
            if (!FoundCode(validCodes, code)) throw std::runtime_error{"invalid code"};

            const auto length = interval.Length();
            auto ref = target.substr(pos, length);
            auto css = query.substr(pos, length);

            const auto NotGap = [](const char c) { return c != '-'; };
            const auto refLen = std::count_if(ref.begin(), ref.end(), NotGap);
            const auto cssLen = std::count_if(ref.begin(), ref.end(), NotGap);

            variant = boost::none;

            switch (code) {
                case 'M':  // fall-through
                case 'N':  // .
                    break;

                case 'R': {
                    assert(css.size() == ref.size());
                    const auto splitupResult = SplitupIUPAC(css);
                    css = splitupResult.readSeq1;
                    variant =
                        Variant{refId,   refPos, refPos + css.size(), ref, splitupResult.readSeq1,
                                refPrev, cssPrev};
                    variant->readSeq2 = splitupResult.readSeq2;
                    variant->frequency1 = splitupResult.freq1;
                    variant->frequency2 = splitupResult.freq2;
                    break;
                }

                case 'I': {
                    const auto splitupResult = SplitupIUPAC(css);
                    css = splitupResult.readSeq1;
                    variant = Variant{refId,   refPos, refPos, "", splitupResult.readSeq1,
                                      refPrev, cssPrev};
                    variant->readSeq2 = splitupResult.readSeq2;
                    variant->frequency1 = splitupResult.freq1;
                    variant->frequency2 = splitupResult.freq2;
                    break;
                }

                case 'D': {
                    variant =
                        Variant{refId, refPos, refPos + ref.size(), ref, "", refPrev, cssPrev};
                    break;
                }

                default:
                    throw std::runtime_error("invalid code");
            }

            if (variant) {
                // HACK ALERT: variants at the first and last position
                // are not handled correctly

                if (!siteCoverage.empty()) {
                    const auto i = std::min(refPos - window.Start(), siteCoverage.size() - 1);
                    variant->coverage = siteCoverage.at(i);
                }
                if (effectiveSiteCoverage && !effectiveSiteCoverage->empty()) {
                    const auto i =
                        std::min(refPos - window.Start(), effectiveSiteCoverage->size() - 1);
                    variant->Annotate("effectiveSiteoverage",
                                      std::to_string(effectiveSiteCoverage->at(i)));
                }
                if (cssQvInWindow && !cssQvInWindow->empty()) {
                    const auto i = std::min(cssPos, cssQvInWindow->size() - 1);
                    variant->confidence = cssQvInWindow->at(i);
                }
                variants.push_back(variant.get());
            }

            // update counters
            refPos += refLen;
            cssPos += cssLen;

            const auto IsGap = [](const char c) { return c == '-'; };
            ref.erase(std::remove_if(ref.begin(), ref.end(), IsGap), ref.end());
            css.erase(std::remove_if(css.begin(), css.end(), IsGap), css.end());

            refPrev = (ref.empty() ? refPrev : ref.back());
            cssPrev = (css.empty() ? cssPrev : css.back());
        }
        return variants;
    }

    static std::string ConstructIUPACFreeConsensus(
        const std::unique_ptr<PacBio::Align::PairwiseAlignment>& ga)
    {
        const auto& target = ga->Target();
        const auto& query = ga->Query();
        assert(target.size() == query.size());

        std::string newCss;

        for (size_t i = 0; i < query.size(); ++i) {
            const auto currentBase = query.at(i);
            if (currentBase != '-') {
                char newBase;
                if (currentBase == 'N' || currentBase == 'n')
                    newBase = currentBase;
                else {
                    const auto targetBase = target.at(i);
                    const std::string iupacCandidates{LookupIUPAC(currentBase)};
                    if (iupacCandidates.find(targetBase) != std::string::npos)
                        newBase = targetBase;
                    else
                        newBase = iupacCandidates.front();
                }
                newCss.push_back(newBase);
            }
        }

        // Be absolutely sure that *really* all ambiguous bases have been removed.
        static constexpr const std::array<char, 12> ambiguousBases{
            {'M', 'm', 'R', 'r', 'W', 'w', 'S', 's', 'Y', 'y', 'K', 'k'}};
        for (const auto c : newCss) {
            if (FoundCode(ambiguousBases, c)) throw std::runtime_error("failed IUPAC resolving");
        }

        return newCss;
    }

    static WindowResult VariantsFromConsensus(
        const ReferenceWindow& window, const std::string& intervalRefSeq, const Consensus& css,
        const std::vector<uint8_t> siteCoverage,
        const boost::optional<std::vector<uint8_t>> effectiveSiteCoverage, const Settings& settings)
    {
        using PacBio::Align::AlignAffine;
        using PacBio::Align::AlignAffineIupac;
        using PairwiseAlignment = PacBio::Align::PairwiseAlignment;

        std::unique_ptr<PairwiseAlignment> ga = nullptr;
        if (settings.polishDiploid)
            ga.reset(AlignAffineIupac(intervalRefSeq, css.sequence));
        else
            ga.reset(AlignAffine(intervalRefSeq, css.sequence));

        std::string newCss;
        if (settings.polishDiploid) {
            newCss = ConstructIUPACFreeConsensus(ga);
            assert(newCss.size() == css.sequence.size());
        }

        return WindowResult{
            Consensus{window, std::move(newCss), std::vector<uint8_t>()},
            VariantsFromAlignment(ga, window, css.confidence, siteCoverage, effectiveSiteCoverage)};
    }

    static std::vector<Variant> FilterVariants(const std::vector<Variant>& variants,
                                               const Settings& settings)
    {
        std::vector<Variant> result;
        for (const auto& v : variants) {
            if (v.coverage.get() >= settings.minCoverage &&
                v.confidence.get() >= settings.minConfidence) {
                result.emplace_back(v);
            }
        }
        return result;
    }

    static void AnnotateVariants(std::vector<Variant>* variants,
                                 const std::vector<PacBio::BAM::BamRecord>& reads)
    {
        for (auto& v : *variants) {
            std::string annotation;
            for (const auto& read : reads) {
                if (!annotation.empty()) {
                    annotation.push_back(',');
                    annotation.push_back(' ');
                }
                annotation.append(read.FullName());
            }
            v.annotations->insert(std::make_pair("rows", annotation));
        }
    }

    static inline void ClipReadsToWindow(std::vector<PacBio::BAM::BamRecord>* reads,
                                         const ReferenceWindow& window)
    {
        const auto winStart = window.Start();
        const auto winEnd = window.End();
        for (auto& read : *reads)
            read.Clip(PacBio::BAM::ClipType::CLIP_TO_REFERENCE, winStart, winEnd);
    }

    static inline void FilterAlignments(std::vector<PacBio::BAM::BamRecord>* reads,
                                        const Settings& settings)
    {
        const auto IsPoaIncompatible = [&](const PacBio::BAM::BamRecord& record) -> bool {
            const auto readLength = record.AlignedEnd() - record.AlignedStart();
            const auto refLength = record.ReferenceEnd() - record.ReferenceStart();
            const auto snr = record.SignalToNoise();
            return (readLength < refLength * settings.readStumpinessThreshold) ||
                   (*std::min_element(snr.begin(), snr.end()) < settings.minHqRegionSnr) ||
                   (record.ReadAccuracy() < settings.minReadScore);
        };

        reads->erase(std::remove_if(reads->begin(), reads->end(), IsPoaIncompatible), reads->end());
    }

    static WindowResult ConsensusAndVariantsForWindow(const Input& input,
                                                      const ReferenceWindow& window,
                                                      const std::string refSeq,
                                                      const Settings& settings)
    {
        using BamRecord = PacBio::BAM::BamRecord;
        using Consensus = PacBio::GenomicConsensus::Consensus;
        using Interval = PacBio::Data::Interval;

        std::vector<Consensus> subconsensi;
        std::vector<Variant> variants;

        // TODO: "fancy chunking"

        std::vector<Interval> intervals{window.interval};
        for (const auto& interval : intervals) {
            // grab interval data
            ReferenceWindow subWindow{window.name, interval};
            const std::string intervalRefSeq = refSeq.substr(interval.Left(), interval.Length());
            auto reads = input.ReadsInWindow(subWindow);
            ClipReadsToWindow(&reads, subWindow);
            FilterAlignments(&reads, settings);

            // if enough POA coverage
            const size_t numSpanning =
                std::count_if(reads.begin(), reads.end(), [&interval](const BamRecord& read) {
                    const auto readStart = static_cast<size_t>(read.ReferenceStart());
                    const auto readEnd = static_cast<size_t>(read.ReferenceEnd());
                    return readStart <= interval.Left() && interval.Right() <= readEnd;
                });

            Consensus css;
            if (numSpanning >= settings.minPoaCoverage) {

                boost::optional<std::vector<uint8_t>> effectiveSiteCoverage;
                if (settings.reportEffectiveCoverage) {
                    std::vector<BamRecord> readsUsed;
                    css =
                        ConsensusForAlignments(window, intervalRefSeq, reads, settings, &readsUsed);
                    effectiveSiteCoverage = CoverageInWindow(subWindow, readsUsed);
                } else
                    css = ConsensusForAlignments(window, intervalRefSeq, reads, settings, nullptr);

                const auto siteCoverage = CoverageInWindow(subWindow, reads);
                const auto tempWindow = VariantsFromConsensus(
                    subWindow, intervalRefSeq, css, siteCoverage, effectiveSiteCoverage, settings);

                auto filteredVariants = FilterVariants(tempWindow.variants, settings);

                if (settings.annotateGFF) AnnotateVariants(&filteredVariants, reads);

                // append filtered variants to final result
                variants.insert(variants.end(), filteredVariants.begin(), filteredVariants.end());

                // The nascent consensus sequence might contain ambiguous bases, these
                // need to be removed as software in the wild cannot deal with such
                // characters and we only use IUPAC for *internal* bookkeeping.
                if (settings.polishDiploid) css.sequence = tempWindow.css.sequence;

                // maybe dump evidence

            } else  // not enough coverage
            {
                css = Consensus::NoCallConsensus(NoCallStyle::LOWERCASE_REFERENCE, subWindow,
                                                 intervalRefSeq);
            }

            // add result to our final consensi
            subconsensi.push_back(css);
        }

        return WindowResult{Consensus::Join(std::move(subconsensi)), std::move(variants)};
    }
};

}  // namespace GenomicConsensus
}  // namespace PacBio
