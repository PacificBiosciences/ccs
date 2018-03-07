// Author: Derek Barnett

#pragma once

#include <cassert>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>

#include <algorithm>
#include <string>
#include <vector>

#include <pbbam/BamRecord.h>
#include <pbcopper/utility/MoveAppend.h>

#include <pacbio/genomicconsensus/experimental/Consensus.h>
#include <pacbio/genomicconsensus/experimental/Input.h>
#include <pacbio/genomicconsensus/experimental/ReferenceWindow.h>
#include <pacbio/genomicconsensus/experimental/Settings.h>
#include <pacbio/genomicconsensus/experimental/Variant.h>
#include <pacbio/genomicconsensus/experimental/WindowResult.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

struct Plurality
{
    static std::string AlignedReference(const PacBio::BAM::Cigar& cigar, const std::string& ref)
    {
        std::string alignedRef;
        alignedRef.reserve(ref.size());

        size_t posInInput = 0;
        for (const auto& op : cigar) {
            const auto len = op.Length();
            if (op.Type() == PacBio::BAM::CigarOperationType::INSERTION)
                alignedRef.append(len, '-');
            else {
                alignedRef.append(ref.substr(posInInput, len));
                posInInput += len;
            }
        }

        return alignedRef;
    }

    using BaseCallsMatrix = std::vector<std::vector<std::string>>;  // mat[read][pos] = bases

    static std::vector<std::string> BaseCallsForAlignment(const std::string& alnRead,
                                                          const std::string& alnRef,
                                                          const size_t windowLength)
    {
        // Idea: scan through the ref, read; for each non-gap character
        // in ref, record all non-gap characters seen in read since
        // last ref gap.

        std::vector<std::string> accum;
        accum.reserve(windowLength);

        std::string readBases;
        for (size_t j = 0; j < alnRef.size(); ++j) {
            const auto refBase = alnRef.at(j);
            const auto readBase = alnRead.at(j);

            if (readBase != '-') readBases.push_back(readBase);

            if (refBase != '-') {
                auto basesForRefPos = (readBases.empty() ? std::string(1, '-') : readBases);
                accum.emplace_back(std::move(basesForRefPos));
                readBases.clear();
            }
        }

        return accum;
    }

    static BaseCallsMatrix TabulateBaseCalls(const Input& input, const ReferenceWindow& window)
    {
        // Go through the reads and build up the structured baseCallsMatrix
        // table, which tabulates the read bases occurring at each reference
        // coordinate in each read.  This code is somewhat tricky, read carefully.

        using ClipType = PacBio::BAM::ClipType;
        using Orientation = PacBio::BAM::Orientation;

        const auto refStart = window.Start();
        const auto refEnd = window.End();

        auto reads = input.ReadsInWindow(window);
        BaseCallsMatrix matrix{reads.size(), std::vector<std::string>(window.Length())};

        for (size_t i = 0; i < reads.size(); ++i) {
            auto& read = reads.at(i);
            read.Clip(ClipType::CLIP_TO_REFERENCE, refStart, refEnd);

            const auto alnRead = read.Sequence(Orientation::GENOMIC, true);
            const auto alnRef = AlignedReference(read.CigarData(), input.ReferenceInWindow(window));
            assert(alnRef.size() == alnRead.size());

            // NOTE (DB): Can't find normalizeHomopolymerGaps() method. Come back later.
            //
            // if realignHomopolymers:
            //     alnRef, alnRead =  normalizeHomopolymerGaps(alnRef, alnRead)

            auto accum = BaseCallsForAlignment(alnRead, alnRef, window.Length());

            const auto s = read.AlignedStart() - refStart;
            const auto e = read.AlignedEnd() - refStart;
            assert((e - s) == accum.size());

            auto& alnBaseCalls = matrix.at(i);
            for (size_t j = 0; j < accum.size(); ++j)
                alnBaseCalls.at(j + s) = std::move(accum.at(j));
        }

        return matrix;
    }

    struct Allele
    {
        std::string bases;
        size_t frequency = 0;

        Allele() = default;
        Allele(std::string bases_, size_t frequency_)
            : bases{std::move(bases)}, frequency{frequency_}
        {
        }
    };

    struct Top2
    {
        Allele firstAllele;
        Allele secondAllele;
        size_t totalCoverage;

        Top2(Allele a1, size_t totalCov) : firstAllele{std::move(a1)}, totalCoverage{totalCov} {}

        Top2(Allele a1, Allele a2, const size_t totalCov)
            : firstAllele{std::move(a1)}, secondAllele{std::move(a2)}, totalCoverage{totalCov}
        {
        }
    };

    static std::vector<Top2> TopAllelesFromMatrix(const BaseCallsMatrix& matrix,
                                                  const size_t windowLength)
    {
        std::vector<Top2> result;
        result.reserve(windowLength);

        std::unordered_map<std::string, size_t> alleleCounts;
        for (size_t pos = 0; pos < windowLength; ++pos) {
            for (size_t row = 0; row < matrix.size(); ++row) {
                assert(matrix.at(row).size() == windowLength);
                auto s = matrix.at(row).at(pos);
                if (!s.empty()) alleleCounts[s]++;
            }

            // calculate totals & save top 2
            Allele top1{};
            Allele top2{};
            size_t totalCoverage = 0;

            if (!alleleCounts.empty()) {
                if (alleleCounts.size() == 1) {
                    const auto a = *alleleCounts.cbegin();
                    top1 = Allele{a.first, a.second};
                    totalCoverage = a.second;
                } else {
                    using MapPair = std::pair<std::string, size_t>;
                    std::vector<MapPair> flattened;
                    flattened.reserve(alleleCounts.size());
                    for (auto& allelesAtPos : alleleCounts) {
                        totalCoverage += allelesAtPos.second;
                        flattened.emplace_back(std::move(allelesAtPos));
                    }
                    assert(flattened.size() > 1);

                    // sort descending on allele count, top2 get placed at beginning
                    std::nth_element(flattened.begin(), flattened.begin() + 1, flattened.end(),
                                     [](const MapPair& lhs, const MapPair& rhs) {
                                         return lhs.second > rhs.second;
                                     });

                    top1 = Allele{std::move(flattened.at(0).first), flattened.at(0).second};
                    top2 = Allele{std::move(flattened.at(1).first), flattened.at(1).second};
                }
            }

            result.emplace_back(top1, top2, totalCoverage);
            alleleCounts.clear();
        }

        return result;
    }

    static std::vector<Top2> TopAllelesForWindow(const Input& input, const ReferenceWindow& window)
    {
        const auto windowLength = window.Length();
        const auto baseCallsMatrix = TabulateBaseCalls(input, window);
        return TopAllelesFromMatrix(baseCallsMatrix, windowLength);
    }

    static std::vector<Variant> VariantsFromRefAndRead(const std::string& refName,
                                                       const size_t refStart, const char refBase,
                                                       const Allele& readAllele,
                                                       const size_t confidence,
                                                       const size_t coverage, const char refPrev,
                                                       const char readPrev)
    {

        // Compute the haploid/heterozygous Variant[s] corresponding to a
        // readSeq aligned against refSeq.
        //
        // Two variant scenario:
        //   REF:   G
        //   READ: AC
        //     => insertion(A), substitution(G->C)
        //
        // Required: refBase != readSeq
        // Returned: List of Variant objects (length one or two)

        assert(std::string(1, refBase) != readAllele.bases);
        assert(!readAllele.bases.empty());

        std::vector<Variant> result;

        const auto& readSeq = readAllele.bases;
        const auto readBefore = readSeq.substr(0, readSeq.size() - 1);

        if (!readBefore.empty()) {
            // Insertion
            Variant v{refName, refStart, refStart, "", readBefore, refPrev, readPrev};
            v.confidence = confidence;
            v.coverage = coverage;
            v.frequency1 = readAllele.frequency;
            result.emplace_back(std::move(v));
        }

        const auto readAt = readSeq.back();
        if (readAt != refBase) {
            Variant v{
                refName, refStart, refStart + 1, std::string(1, refBase), std::string(1, readAt),
                refPrev, readPrev};
            v.confidence = confidence;
            v.coverage = coverage;
            v.frequency1 = readAllele.frequency;
            result.emplace_back(std::move(v));
        }

        return result;
    }

    static std::vector<Variant> VariantsFromRefAndReads(
        const std::string& refName, const size_t refStart, const char refBase,
        const Allele& cssAllele, const Allele& altAllele, const size_t confidence,
        const size_t coverage, const char refPrev, const char readPrev)
    {
        //
        // Heterozygous extension of the above
        //
        assert(std::string(1, refBase) != cssAllele.bases);
        assert(std::string(1, refBase) != altAllele.bases);
        assert(!cssAllele.bases.empty());
        assert(!altAllele.bases.empty());

        std::vector<Variant> result;

        const auto& cssSeq = cssAllele.bases;
        const auto& altSeq = altAllele.bases;
        const auto cssBefore = cssSeq.substr(0, cssSeq.size() - 1);
        const auto altBefore = altSeq.substr(0, altSeq.size() - 1);

        if (!cssBefore.empty() || !altBefore.empty()) {
            // Insertion
            Variant v{refName, refStart, refStart, "", cssBefore, refPrev, readPrev};
            v.readSeq2 = altBefore;
            v.confidence = confidence;
            v.coverage = coverage;
            v.frequency1 = cssAllele.frequency;
            v.frequency2 = altAllele.frequency;
            result.emplace_back(std::move(v));
        }

        const char cssAt = cssSeq.back();
        const char altAt = altSeq.back();
        if (cssAt != refBase || altAt != refBase) {
            Variant v{
                refName, refStart, refStart + 1, std::string(1, refBase), std::string(1, cssAt),
                refPrev, readPrev};
            v.readSeq2 = std::string(1, altAt);
            v.confidence = confidence;
            v.coverage = coverage;
            v.frequency1 = cssAllele.frequency;
            v.frequency2 = altAllele.frequency;
            result.emplace_back(std::move(v));
        }

        return result;
    }

    static bool IsAllUpper(const std::string& seq)
    {
        if (seq.empty()) return false;
        return !std::any_of(seq.cbegin(), seq.cend(), [](const char c) { return ::islower(c); });
    }

    static std::vector<Variant> ComputeVariants(
        const ReferenceWindow& window, const std::string& refSeq,
        const std::vector<size_t>& effectiveCoverage, const std::vector<Allele>& consensusAlleles,
        const std::vector<uint8_t>& consensusConfidences,
        const std::vector<Allele>& alternateAlleles,          // DIPLOID ONLY
        const std::vector<uint8_t>& heterozygousConfidences,  // "
        const Settings& settings)
    {
        std::vector<Variant> result;

        const auto& refName = window.name;
        const auto refStart = window.Start();
        const auto windowSize = window.Length();

        assert(refSeq.size() == windowSize);
        assert(consensusAlleles.size() == windowSize);

        if (settings.diploid) assert(alternateAlleles.size() == windowSize);

        char refPrev = 'N';
        char cssPrev = 'N';
        for (size_t j = 0; j < windowSize; ++j) {
            const auto cov = effectiveCoverage.at(j);
            if (cov < settings.minCoverage) continue;

            const auto refPos = refStart + j;
            const auto refBase = refSeq.at(j);
            const auto cssAllele = consensusAlleles.at(j);
            const auto conf = consensusConfidences.at(j);

            const size_t hetConf = heterozygousConfidences.at(j);
            if (settings.diploid && hetConf > conf) {
                if (hetConf >= settings.minConfidence && refBase != 'N') {
                    const auto& altAllele = alternateAlleles.at(j);
                    auto vars = VariantsFromRefAndReads(refName, refPos, refBase, cssAllele,
                                                        altAllele, hetConf, cov, refPrev, cssPrev);
                    PacBio::Utility::MoveAppend(std::move(vars), result);
                }
            } else {
                //
                // Haploid variant[s]?
                //
                if ((conf >= settings.minConfidence) && (refBase != 'N') &&
                    (cssAllele.bases != "N") &&
                    (cssAllele.bases.empty() || IsAllUpper(cssAllele.bases)) &&
                    (std::string(1, refBase) != cssAllele.bases)) {
                    auto vars = VariantsFromRefAndRead(refName, refPos, refBase, cssAllele, conf,
                                                       cov, refPrev, cssPrev);
                    PacBio::Utility::MoveAppend(std::move(vars), result);
                }
            }

            // if we have ref or css bases, update the anchors
            refPrev = (refBase == '\0' ? refPrev : refBase);
            cssPrev = (cssAllele.bases.empty() ? cssPrev : cssAllele.bases.back());
        }

        // filter, sort, & return result
        if (settings.diploid) {
            result.erase(
                std::remove_if(result.begin(), result.end(),
                               [](const Variant& v) {
                                   const auto altSize =
                                       (v.readSeq2 == boost::none ? 0 : v.readSeq2->size());
                                   return v.IsHeterozygous() && v.readSeq1.size() != altSize;
                               }),
                result.end());
        }

        std::sort(result.begin(), result.end());
        return result;
    }

    struct SiteConfidences
    {
        uint8_t consensusConfidence;
        uint8_t heterozygousConfidence;
    };

    static SiteConfidences PosteriorConfidences(size_t depth, size_t cssFreq, size_t altFreq,
                                                const bool diploid = false, const double cap = 40)
    {
        // FIXME (DB): math here is off when (diploid == true)

        //
        // Return crude approximations to the posterior probabilities of the
        // genotypes s_1 and s_1/s_2, where s_1 and s_2 are the observed
        // consensus and alternate allele.  The assumption here is that the
        // probability of the genotype being anything other that s_1, s_2, or
        // s_1/s_2 is vanishingly small.  Not really a very good assumption,
        // but plurality is not our real algorithm anyway.
        //

        static constexpr const double EPS = 0.05;
        static const double LOGEPS = std::log(EPS);
        static const double LOG_O_M_EPS = std::log(1 - EPS);
        static const double LOG_O_M_EPS_2 = std::log((1 - EPS) / 2);

        ++cssFreq;
        ++altFreq;
        depth += 2;

        const double cssLL = (cssFreq * LOG_O_M_EPS) + ((depth - cssFreq) * LOGEPS);
        const double altLL = (altFreq * LOG_O_M_EPS) + ((depth - altFreq) * LOGEPS);
        const double cssL = std::exp(cssLL);
        const double altL = std::exp(altLL);

        double total = 0.0;
        double hetConf = 0.0;
        if (diploid) {
            const auto hetLL =
                ((cssFreq + altFreq) * LOG_O_M_EPS_2) + ((depth - cssFreq - altFreq) * LOGEPS);
            const auto hetL = std::exp(hetLL);
            total = cssL + altL + hetL;

            const auto hetProb = hetL / total;
            hetConf = ((hetProb < 1.0) ? (-10 * std::log10(1.0 - hetProb)) : cap);
        } else {
            total = cssL + altL;
            hetConf = 0.0;
        }

        const double cssProb = cssL / total;
        const double cssConf = ((cssProb < 1.0) ? (-10 * std::log10(1.0 - cssProb)) : cap);

        return SiteConfidences{static_cast<uint8_t>(std::min(cap, cssConf)),
                               static_cast<uint8_t>(std::min(cap, hetConf))};
    }

    static WindowResult ConsensusAndVariantsForWindow(const Input& input,
                                                      const ReferenceWindow& window,
                                                      const std::string refSeq,
                                                      const Settings& settings)
    {
        //
        // Compute (Consensus, [Variant]) for this window, using the given
        // `alns`, by applying a straightforward column-oriented consensus
        // calling algorithm.
        //
        // If the consensus cannot be called for a base, "N" will be placed
        // in the consensus sequence for that position.
        //
        // If `realignHomopolymers` is True, alignment gaps will be shuffled
        // in homopolymer regions in an attempt to maximize variant detection
        // sensitivity (not yet implemented, and may never be).
        //

        //
        // Build up these arrays in reference coordinates.
        //

        std::vector<Allele> consensusAlleles;
        std::vector<uint8_t> consensusConfidences;
        std::vector<size_t> effectiveCoverage;
        std::vector<Allele> alternateAlleles;          // DIPLOID ONLY
        std::vector<uint8_t> heterozygousConfidences;  // ""

        const auto noCallConsensus =
            Consensus::NoCallConsensus(NoCallStyle::NO_CALL, window, refSeq);
        const auto topAllelesForWindow = TopAllelesForWindow(input, window);

        for (size_t i = 0; i < window.Length(); ++i) {
            const auto topAlleles = topAllelesForWindow.at(i);

            Allele siteConsensusAllele;
            const auto siteEffectiveCoverage = topAlleles.totalCoverage;
            if ((siteEffectiveCoverage == 0) || (siteEffectiveCoverage > settings.minCoverage)) {
                siteConsensusAllele =
                    Allele{std::string(1, noCallConsensus.sequence.at(i)), siteEffectiveCoverage};
            } else
                siteConsensusAllele = topAlleles.firstAllele;

            // Replace explicit gaps with empty string
            if (siteConsensusAllele.bases == "-") siteConsensusAllele.bases.clear();

            consensusAlleles.push_back(siteConsensusAllele);
            effectiveCoverage.push_back(siteEffectiveCoverage);

            Allele siteAlternateAllele;
            if (settings.diploid) {
                if (topAlleles.secondAllele.frequency > 0)
                    siteAlternateAllele = topAlleles.secondAllele;
                if (siteAlternateAllele.bases == "-") siteAlternateAllele.bases.clear();
                alternateAlleles.push_back(siteAlternateAllele);
            }

            const auto siteConfidences =
                PosteriorConfidences(siteEffectiveCoverage, siteConsensusAllele.frequency,
                                     siteAlternateAllele.frequency, settings.diploid);

            consensusConfidences.push_back(siteConfidences.consensusConfidence);
            if (settings.diploid)
                heterozygousConfidences.push_back(siteConfidences.heterozygousConfidence);
        }

        //
        // Derive variants from reference-coordinates consensus
        //
        auto variants = ComputeVariants(window, refSeq, effectiveCoverage, consensusAlleles,
                                        consensusConfidences, alternateAlleles,
                                        heterozygousConfidences, settings);

        //
        // Now we need to put everything in consensus coordinates
        //
        std::string consensusSequence;
        std::vector<uint8_t> consensusConfidenceOut;
        for (size_t i = 0; i < consensusAlleles.size(); ++i) {
            auto& allele = consensusAlleles.at(i);
            const auto alleleLen = static_cast<uint8_t>(allele.bases.size());

            consensusSequence.append(allele.bases);
            PacBio::Utility::MoveAppend(std::vector<uint8_t>{alleleLen, consensusConfidences.at(i)},
                                        consensusConfidenceOut);
        }

        return WindowResult{
            Consensus{window, std::move(consensusSequence), std::move(consensusConfidenceOut)},
            std::move(variants)};
    }
};

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
