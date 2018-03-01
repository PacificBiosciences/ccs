// Author: Derek Barnett

#include <pacbio/align/BandedChainAlignment.h>

#include <cassert>
#include <cfloat>

#include <algorithm>
#include <array>
#include <utility>
#include <vector>

#include <pacbio/align/internal/BCAlignBlocks.h>
#include <pacbio/align/internal/BCAlignImpl.h>

namespace PacBio {
namespace Align {
namespace Internal {

static inline float max4(const float a, const float b, const float c, const float d)
{
    return std::max(std::max(a, b), std::max(c, d));
}

static inline float score(const char t, const char q, const BandedChainAlignConfig& config)
{
    return (t == q ? config.matchScore_ : config.mismatchPenalty_);
}

static inline void addCigarOp(PacBio::Data::Cigar* cigar, const PacBio::Data::CigarOperationType op)
{
    if (cigar->empty())
        cigar->emplace_back(op, 1);
    else {
        auto& lastOp = cigar->back();
        if (lastOp.Type() == op)
            lastOp.Length(lastOp.Length() + 1);  // extend current CIGAR op length
        else
            cigar->emplace_back(op, 1);
    }
}

// ------------------------
// BandedGlobalAlignBlock
// ------------------------

PacBio::Data::Cigar BandedGlobalAlignBlock::Align(const char* target, const char* query,
                                                  PacBio::Align::Seed seed)
{
    using PacBio::Data::Cigar;
    using PacBio::Data::CigarOperation;
    using PacBio::Data::CigarOperationType;

    Cigar cigar;

    // ensure horizontal sequence length is >= vertical
    // (simplifies band calculations)
    const size_t qLen = seed.EndPositionV() - seed.BeginPositionV();
    const size_t tLen = seed.EndPositionH() - seed.BeginPositionH();
    if (qLen == 0) {
        cigar.emplace_back(CigarOperationType::DELETION, tLen);
        return cigar;
    } else if (tLen == 0) {
        cigar.emplace_back(CigarOperationType::INSERTION, tLen);
        return cigar;
    }

    const bool seqsFlipped = (qLen > tLen);
    const char* seq1 =
        (seqsFlipped ? (target + seed.BeginPositionH()) : (query + seed.BeginPositionV()));
    const char* seq2 =
        (seqsFlipped ? (query + seed.BeginPositionV()) : (target + seed.BeginPositionH()));
    const size_t seq1Len = (seqsFlipped ? tLen : qLen);
    const size_t seq2Len = (seqsFlipped ? qLen : tLen);

    // Initialize space & scores
    Init(seq2Len, seq1Len);

    // for each row
    for (size_t i = 1; i <= seq1Len; ++i) {

        // foreach column within band
        const auto& e = lookup_.at(i);
        for (size_t j = e.jBegin_; j <= e.jEnd_; ++j) {

            if (j == 0) continue;

            const auto currentIdx = IndexFor(i, j);
            const auto diagIdx = IndexFor(i - 1, j - 1);
            const auto upIdx = IndexFor(i - 1, j);
            const auto leftIdx = IndexFor(i, j - 1);
            const bool upAllowed = (upIdx != std::string::npos);
            const bool leftAllowed = (leftIdx != std::string::npos);

            const auto s = score(seq2[j - 1], seq1[i - 1], config_);

            matchScores_.at(currentIdx) =
                (std::max(matchScores_.at(diagIdx), gapScores_.at(diagIdx)) + s);

            gapScores_.at(currentIdx) =
                max4((leftAllowed ? matchScores_.at(leftIdx) + config_.gapOpenPenalty_ : -FLT_MAX),
                     (leftAllowed ? gapScores_.at(leftIdx) + config_.gapExtendPenalty_ : -FLT_MAX),
                     (upAllowed ? matchScores_.at(upIdx) + config_.gapOpenPenalty_ : -FLT_MAX),
                     (upAllowed ? gapScores_.at(upIdx) + config_.gapExtendPenalty_ : -FLT_MAX));
        }
    }

    // Traceback
    const size_t MATCH_MATRIX = 1;
    const size_t GAP_MATRIX = 2;

    // find traceback start
    const auto btStart = BacktraceStart(seq2Len, seq1Len);
    size_t i = btStart.first;
    size_t j = btStart.second;
    const size_t backtraceStartIdx = IndexFor(btStart.first, btStart.second);

    size_t mat =
        (matchScores_.at(backtraceStartIdx) >= gapScores_.at(backtraceStartIdx) ? MATCH_MATRIX
                                                                                : GAP_MATRIX);
    size_t iPrev, jPrev, matPrev;

    // if not beginning at bottom right)
    if (i < seq1Len) {
        const auto op =
            (seqsFlipped ? CigarOperationType::DELETION : CigarOperationType::INSERTION);
        for (size_t k = seq1Len - i; k > 0; --k)
            addCigarOp(&cigar, op);
    } else if (j < seq2Len) {
        const auto op =
            (seqsFlipped ? CigarOperationType::INSERTION : CigarOperationType::DELETION);
        for (size_t k = seq2Len - j; k > 0; --k)
            addCigarOp(&cigar, op);
    }

    while (i > 0 || j > 0) {

        if (mat == MATCH_MATRIX) {

            const auto diagIdx = IndexFor(i - 1, j - 1);
            matPrev =
                (matchScores_.at(diagIdx) >= gapScores_.at(diagIdx) ? MATCH_MATRIX : GAP_MATRIX);
            iPrev = i - 1;
            jPrev = j - 1;
            const auto op = (seq1[iPrev] == seq2[jPrev] ? CigarOperationType::SEQUENCE_MATCH
                                                        : CigarOperationType::SEQUENCE_MISMATCH);
            addCigarOp(&cigar, op);

        } else {
            assert(mat == GAP_MATRIX);

            const auto upIdx = IndexFor(i - 1, j);
            const auto leftIdx = IndexFor(i, j - 1);
            const auto upAllowed = (upIdx != std::string::npos);
            const auto leftAllowed = (leftIdx != std::string::npos);

            const std::array<float, 4> s{
                {(j > 0 && leftAllowed ? matchScores_.at(leftIdx) + config_.gapOpenPenalty_
                                       : -FLT_MAX),
                 (j > 0 && leftAllowed ? gapScores_.at(leftIdx) + config_.gapExtendPenalty_
                                       : -FLT_MAX),
                 (i > 0 && upAllowed ? matchScores_.at(upIdx) + config_.gapOpenPenalty_ : -FLT_MAX),
                 (i > 0 && upAllowed ? gapScores_.at(upIdx) + config_.gapExtendPenalty_
                                     : -FLT_MAX)}};
            const auto argMax = std::distance(s.cbegin(), std::max_element(s.cbegin(), s.cend()));

            matPrev = ((argMax == 0 || argMax == 2) ? MATCH_MATRIX : GAP_MATRIX);
            if (argMax == 0 || argMax == 1) {

                iPrev = i;
                jPrev = j - 1;
                const auto op =
                    (seqsFlipped ? CigarOperationType::INSERTION : CigarOperationType::DELETION);
                addCigarOp(&cigar, op);
            } else {

                iPrev = i - 1;
                jPrev = j;
                const auto op =
                    (seqsFlipped ? CigarOperationType::DELETION : CigarOperationType::INSERTION);
                addCigarOp(&cigar, op);
            }
        }

        // step back one
        i = iPrev;
        j = jPrev;
        mat = matPrev;

        // if at edge, force gap moves from here on
        if (i == 0 || j == 0) mat = GAP_MATRIX;
    }

    // reverse CIGAR & return
    std::reverse(cigar.begin(), cigar.end());
    return cigar;
}

std::pair<size_t, size_t> BandedGlobalAlignBlock::BacktraceStart(const size_t tLen,
                                                                 const size_t qLen) const
{
    // NOTE: Finding backtrace start this way allows us to not penalize end-gaps.
    const size_t maxIndex = std::min(qLen, tLen);

    // find max score in last column
    std::pair<size_t, size_t> maxCellRight{maxIndex, maxIndex};
    float maxScoreRight = -FLT_MAX;
    {
        for (size_t i = 1; i <= maxIndex; ++i) {
            const auto& e = lookup_.at(i);
            const size_t lastColumn = e.jEnd_;
            const auto idx = IndexFor(i, lastColumn);
            if (matchScores_[idx] > maxScoreRight) {
                maxScoreRight = matchScores_[idx];
                maxCellRight = std::make_pair(i, lastColumn);
            }
        }
    }

    // find max score in last row
    std::pair<size_t, size_t> maxCellBottom{maxIndex, maxIndex};
    float maxScoreBottom = -FLT_MAX;
    {
        const size_t lastRow = maxIndex;
        const auto& lookupElement = lookup_.at(lastRow);
        for (size_t j = lookupElement.jBegin_; j < lookupElement.jEnd_; ++j) {
            const auto idx = IndexFor(lastRow, j);
            if (matchScores_[idx] > maxScoreBottom) {
                maxScoreBottom = matchScores_[idx];
                maxCellBottom = std::make_pair(lastRow, j);
            }
        }
    }

    return (maxScoreBottom > maxScoreRight ? maxCellBottom : maxCellRight);
}

size_t BandedGlobalAlignBlock::IndexFor(const size_t i, const size_t j) const
{
    // if in matrix
    if (i != std::string::npos && j != std::string::npos) {
        const auto& e = lookup_.at(i);
        // return array index for (i,j) if in-band
        if (j >= e.jBegin_ && j <= e.jEnd_) return e.arrayStart_ + (j - e.jBegin_);
    }

    // (i,j) either out of matrix bounds, or out-of-band
    return std::string::npos;
}

void BandedGlobalAlignBlock::Init(const size_t tLen, const size_t qLen)
{
    const auto numElements = InitLookup(tLen, qLen);
    InitScores(tLen, qLen, numElements);
}

size_t BandedGlobalAlignBlock::InitLookup(const size_t tLen, const size_t qLen)
{
    // ensure space
    lookup_.clear();
    lookup_.reserve(qLen + 1);

    const auto t = static_cast<int64_t>(tLen);
    const auto q = static_cast<int64_t>(qLen);
    const auto k = static_cast<int64_t>(config_.bandExtend_);
    assert(t >= q);

    size_t arrayStart = 0;
    int64_t jBegin = 0;
    int64_t jEnd = 0;
    for (int64_t i = 0; i < q + 1; ++i) {

        // jBegin
        jBegin = i - k;
        if (jBegin < 0) jBegin = 0;

        // jEnd
        jEnd = i + k;
        if (jEnd > t) jEnd = t;

        // sanity checks
        assert(jBegin >= 0);
        assert(jEnd >= 0);
        assert(jEnd >= jBegin);

        // store lookup values
        lookup_.emplace_back(arrayStart, static_cast<size_t>(jBegin), static_cast<size_t>(jEnd));

        // update arrayStart for next row (or 'numElements' on exit)
        const auto rowElements = static_cast<size_t>(jEnd - jBegin) + 1;
        arrayStart += (rowElements);
    }

    return arrayStart;
}

void BandedGlobalAlignBlock::InitScores(const size_t tLen, const size_t qLen, const size_t n)
{
    matchScores_.resize(n);
    gapScores_.resize(n);

    matchScores_[0] = 0;
    gapScores_[0] = -FLT_MAX;

    const auto maxQ = std::min(qLen, config_.bandExtend_);
    const auto maxT = std::min(tLen, config_.bandExtend_);

    for (size_t i = 1; i <= maxQ; ++i) {
        const auto idx = IndexFor(i, 0);
        matchScores_[idx] = -FLT_MAX;
        gapScores_[idx] = config_.gapOpenPenalty_ + (i - 1) * config_.gapExtendPenalty_;
    }

    for (size_t j = 1; j <= maxT; ++j) {
        const auto idx = IndexFor(0, j);
        matchScores_[idx] = -FLT_MAX;
        gapScores_[idx] = config_.gapOpenPenalty_ + (j - 1) * config_.gapExtendPenalty_;
    }
}

// --------------------------
// StandardGlobalAlignBlock
// --------------------------

PacBio::Data::Cigar StandardGlobalAlignBlock::Align(const char* target, const size_t tLen,
                                                    const char* query, const size_t qLen)
{
    using PacBio::Data::Cigar;
    using PacBio::Data::CigarOperation;
    using PacBio::Data::CigarOperationType;

    // Initialize space & scores
    Init(tLen, qLen);

    // Main loop
    for (size_t i = 1; i <= qLen; ++i) {
        for (size_t j = 1; j <= tLen; ++j) {
            const auto s = score(target[j - 1], query[i - 1], config_);
            matchScores_[i][j] =
                (std::max(matchScores_[i - 1][j - 1], gapScores_[i - 1][j - 1]) + s);
            gapScores_[i][j] = max4(matchScores_[i][j - 1] + config_.gapOpenPenalty_,
                                    gapScores_[i][j - 1] + config_.gapExtendPenalty_,
                                    matchScores_[i - 1][j] + config_.gapOpenPenalty_,
                                    gapScores_[i - 1][j] + config_.gapExtendPenalty_);
        }
    }

    // Traceback
    const size_t MATCH_MATRIX = 1;
    const size_t GAP_MATRIX = 2;

    // find traceback start
    const auto btStart = BacktraceStart(tLen, qLen);
    size_t i = btStart.first;
    size_t j = btStart.second;
    size_t mat =
        (matchScores_[btStart.first][btStart.second] >= gapScores_[btStart.first][btStart.second]
             ? MATCH_MATRIX
             : GAP_MATRIX);
    size_t iPrev, jPrev, matPrev;
    Cigar cigar;

    // if not beginning at bottom right, add corresponding indel
    if (i < qLen) {
        for (size_t k = qLen - i; k > 0; --k)
            addCigarOp(&cigar, CigarOperationType::INSERTION);
    } else if (j < tLen) {
        for (size_t k = tLen - j; k > 0; --k)
            addCigarOp(&cigar, CigarOperationType::DELETION);
    }

    // traceback remaining sequence
    while (i > 0 || j > 0) {

        if (mat == MATCH_MATRIX) {
            matPrev = (matchScores_[i - 1][j - 1] >= gapScores_[i - 1][j - 1] ? MATCH_MATRIX
                                                                              : GAP_MATRIX);
            iPrev = i - 1;
            jPrev = j - 1;
            const auto op = (query[iPrev] == target[jPrev] ? CigarOperationType::SEQUENCE_MATCH
                                                           : CigarOperationType::SEQUENCE_MISMATCH);
            addCigarOp(&cigar, op);

        } else {
            assert(mat == GAP_MATRIX);

            const std::array<float, 4> s{
                {(j > 0 ? matchScores_[i][j - 1] + config_.gapOpenPenalty_ : -FLT_MAX),
                 (j > 0 ? gapScores_[i][j - 1] + config_.gapExtendPenalty_ : -FLT_MAX),
                 (i > 0 ? matchScores_[i - 1][j] + config_.gapOpenPenalty_ : -FLT_MAX),
                 (i > 0 ? gapScores_[i - 1][j] + config_.gapExtendPenalty_ : -FLT_MAX)}};
            const auto argMax = std::distance(s.cbegin(), std::max_element(s.cbegin(), s.cend()));

            matPrev = ((argMax == 0 || argMax == 2) ? MATCH_MATRIX : GAP_MATRIX);
            if (argMax == 0 || argMax == 1) {
                iPrev = i;
                jPrev = j - 1;
                addCigarOp(&cigar, CigarOperationType::DELETION);
            } else {
                iPrev = i - 1;
                jPrev = j;
                addCigarOp(&cigar, CigarOperationType::INSERTION);
            }
        }

        // step back one
        i = iPrev;
        j = jPrev;
        mat = matPrev;
    }

    // reverse CIGAR & return
    std::reverse(cigar.begin(), cigar.end());
    return cigar;
}

std::pair<size_t, size_t> StandardGlobalAlignBlock::BacktraceStart(const size_t tLen,
                                                                   const size_t qLen) const
{
    // NOTE: Finding backtrace start this way allows us to not penalize end-gaps.

    // find max score in last column
    std::pair<size_t, size_t> maxCellRight{qLen, tLen};
    float maxScoreRight = -FLT_MAX;
    const size_t lastColumn = tLen;
    for (size_t i = 1; i <= qLen; ++i) {
        if (matchScores_[i][lastColumn] > maxScoreRight) {
            maxScoreRight = matchScores_[i][lastColumn];
            maxCellRight = std::make_pair(i, lastColumn);
        }
    }

    // find max score in last row
    std::pair<size_t, size_t> maxCellBottom{qLen, tLen};
    float maxScoreBottom = -FLT_MAX;
    const size_t lastRow = qLen;
    for (size_t j = 1; j <= tLen; ++j) {
        if (matchScores_[lastRow][j] > maxScoreBottom) {
            maxScoreBottom = matchScores_[lastRow][j];
            maxCellBottom = std::make_pair(lastRow, j);
        }
    }

    return (maxScoreBottom > maxScoreRight ? maxCellBottom : maxCellRight);
}

void StandardGlobalAlignBlock::Init(const size_t tLen, const size_t qLen)
{
    // ensure space
    matchScores_.resize(qLen + 1);
    gapScores_.resize(qLen + 1);

    assert(matchScores_.size() == qLen + 1);
    assert(gapScores_.size() == qLen + 1);

    for (size_t i = 0; i < qLen + 1; ++i) {
        matchScores_[i].resize(tLen + 1);
        gapScores_[i].resize(tLen + 1);
    }

    // fill out initial scores
    matchScores_[0][0] = 0;
    gapScores_[0][0] = -FLT_MAX;
    for (size_t i = 1; i <= qLen; ++i) {
        matchScores_[i][0] = -FLT_MAX;
        gapScores_[i][0] = config_.gapOpenPenalty_ + (i - 1) * config_.gapExtendPenalty_;
    }
    for (size_t j = 1; j <= tLen; ++j) {
        matchScores_[0][j] = -FLT_MAX;
        gapScores_[0][j] = config_.gapOpenPenalty_ + (j - 1) * config_.gapExtendPenalty_;
    }
}

// ------------------------
// BandedChainAlignerImpl
// ------------------------

using SeedVector = std::vector<PacBio::Align::Seed>;
using SeedsIter = SeedVector::iterator;
using SeedsConstIter = SeedVector::const_iterator;

static SeedsConstIter FirstAnchorSeed(const SeedVector& seeds, const size_t band)
{
    return std::find_if(seeds.cbegin(), seeds.cend(), [band](const PacBio::Align::Seed& seed) {
        return seed.BeginPositionH() >= band && seed.BeginPositionV() >= band;
    });
}

static SeedsConstIter LastAnchorSeed(const SeedVector& seeds, const size_t tLen, const size_t qLen,
                                     const size_t band)
{
    return std::find_if(seeds.crbegin(), seeds.crend(),
                        [tLen, qLen, band](const PacBio::Align::Seed& seed) {
                            return seed.EndPositionH() + band < tLen &&
                                   seed.EndPositionV() + band < qLen;
                        })
        .base();  // convert reverse_iterator -> iterator
}

BandedChainAlignerImpl::BandedChainAlignerImpl(const BandedChainAlignConfig& config)
    : config_(config), gapBlock_(config), seedBlock_(config), gapBlockBeginH_(0), gapBlockBeginV_(0)
{
}

BandedChainAlignment BandedChainAlignerImpl::Align(const char* target, const size_t targetLen,
                                                   const char* query, const size_t queryLen,
                                                   const std::vector<PacBio::Align::Seed>& seeds)
{
    // return empty alignment on empty seeds
    if (seeds.empty()) return BandedChainAlignment{};

    // reset state & store input sequence info
    Initialize(target, targetLen, query, queryLen);

    // step through merged seeds (all overlaps collapsed)
    //   1 - align gap region before current seed, and then
    //   2 - align current seed
    const auto mergedSeeds = MergeSeeds(seeds);
    const auto band = config_.bandExtend_;
    auto it = FirstAnchorSeed(mergedSeeds, band);
    const auto end = LastAnchorSeed(mergedSeeds, targetLen, queryLen, band);
    for (; it != end; ++it) {
        AlignGapBlock(*it);
        AlignSeedBlock(*it);
    }

    // finally align last gap region after last seed & return result
    AlignLastGapBlock();
    return Result();
}

void BandedChainAlignerImpl::AlignGapBlock(const size_t hLength, const size_t vLength)
{
    // do 'standard' DP align
    auto cigar = gapBlock_.Align(sequences_.target + gapBlockBeginH_, hLength,
                                 sequences_.query + gapBlockBeginV_, vLength);

    // incorporate alignment into total result
    StitchCigars(&globalCigar_, std::move(cigar));
}

void BandedChainAlignerImpl::AlignGapBlock(const PacBio::Align::Seed& nextSeed)
{
    const size_t hLength = nextSeed.BeginPositionH() - gapBlockBeginH_;
    const size_t vLength = nextSeed.BeginPositionV() - gapBlockBeginV_;
    AlignGapBlock(hLength, vLength);
}

void BandedChainAlignerImpl::AlignLastGapBlock()
{
    const size_t hLength = sequences_.targetLen - gapBlockBeginH_;
    const size_t vLength = sequences_.queryLen - gapBlockBeginV_;
    AlignGapBlock(hLength, vLength);
}

void BandedChainAlignerImpl::AlignSeedBlock(const PacBio::Align::Seed& seed)
{
    // do seed-guided, banded align
    auto cigar = seedBlock_.Align(sequences_.target, sequences_.query, seed);

    // incorporate alignment into total result
    StitchCigars(&globalCigar_, std::move(cigar));

    // see if we ended with an indel, if so remove that and try re-aligning that
    // portion in the next alignment phase
    using PacBio::Data::CigarOperationType;
    size_t hOffset = 0;
    size_t vOffset = 0;
    const auto& lastOp = globalCigar_.back();
    if (lastOp.Type() == CigarOperationType::DELETION) {
        hOffset = lastOp.Length();
        globalCigar_.pop_back();
    } else if (lastOp.Type() == CigarOperationType::INSERTION) {
        vOffset = lastOp.Length();
        globalCigar_.pop_back();
    }

    // update offsets for next alignment block
    gapBlockBeginH_ = seed.EndPositionH() - hOffset;
    gapBlockBeginV_ = seed.EndPositionV() - vOffset;
}

void BandedChainAlignerImpl::Initialize(const char* target, const size_t targetLen,
                                        const char* query, const size_t queryLen)
{
    globalCigar_.clear();
    globalScore_ = std::numeric_limits<int64_t>::min();
    gapBlockBeginH_ = 0;
    gapBlockBeginV_ = 0;

    sequences_ = Sequences{target, targetLen, query, queryLen};
}

std::vector<PacBio::Align::Seed> BandedChainAlignerImpl::MergeSeeds(
    const std::vector<PacBio::Align::Seed>& seeds)
{
    // no merging needed on empty or single-element containers
    if (seeds.size() <= 1) return seeds;

    // push first seed into output container
    std::vector<PacBio::Align::Seed> mergedSeeds;
    mergedSeeds.reserve(seeds.size());
    mergedSeeds.push_back(seeds.front());
    auto currentSeed = mergedSeeds.begin();

    // iterate over remaining seeds
    auto inputIter = seeds.cbegin() + 1;
    auto inputEnd = seeds.cend();
    for (; inputIter != inputEnd; ++inputIter) {

        // if input seed starts after current seed
        if (inputIter->BeginPositionH() > currentSeed->EndPositionH() &&
            inputIter->BeginPositionV() > currentSeed->EndPositionV()) {
            mergedSeeds.push_back(*inputIter);
            currentSeed = mergedSeeds.end() - 1;
        }

        // else contiguous or overlapping
        else {
            currentSeed->BeginPositionH(
                std::min(inputIter->BeginPositionH(), currentSeed->BeginPositionH()));
            currentSeed->EndPositionH(
                std::max(inputIter->EndPositionH(), currentSeed->EndPositionH()));
            currentSeed->BeginPositionV(
                std::min(inputIter->BeginPositionV(), currentSeed->BeginPositionV()));
            currentSeed->EndPositionV(
                std::max(inputIter->EndPositionV(), currentSeed->EndPositionV()));
        }
    }
    return mergedSeeds;
}

BandedChainAlignment BandedChainAlignerImpl::Result()
{
    return BandedChainAlignment{config_,          sequences_.target,   sequences_.targetLen,
                                sequences_.query, sequences_.queryLen, globalCigar_};
}

void BandedChainAlignerImpl::StitchCigars(PacBio::Data::Cigar* global, PacBio::Data::Cigar&& local)
{
    assert(global);

    // quick checks if either CIGAR empty
    if (local.empty()) return;
    if (global->empty()) {
        *global = std::move(local);
        return;
    }

    global->reserve(global->size() + local.size());

    // see if we can merge first local CIGAR op into last global op
    auto& firstLocalOp = local.front();
    auto& lastGlobalOp = global->back();
    size_t i = 0;
    if (firstLocalOp.Type() == lastGlobalOp.Type()) {
        lastGlobalOp.Length(lastGlobalOp.Length() + firstLocalOp.Length());
        ++i;
    }

    // append remaining local ops to global CIGAR
    for (; i < local.size(); ++i)
        global->emplace_back(std::move(local.at(i)));
}

}  // namespace internal

// ------------------------
// BandedChainAlignment
// ------------------------

BandedChainAlignment::BandedChainAlignment(const BandedChainAlignConfig& config, std::string target,
                                           std::string query, PacBio::Data::Cigar cigar)
    : config_(config)
    , target_{std::move(target)}
    , query_{std::move(query)}
    , cigar_{std::move(cigar)}
{
    using PacBio::Data::CigarOperationType;

    alignedTarget_.reserve(target_.size());
    alignedQuery_.reserve(query_.size());

    size_t tPos = 0;
    size_t qPos = 0;
    for (const auto& op : cigar_) {
        const auto type = op.Type();
        for (size_t i = 0; i < op.Length(); ++i) {
            switch (type) {
                case CigarOperationType::SEQUENCE_MATCH:
                case CigarOperationType::SEQUENCE_MISMATCH:
                case CigarOperationType::ALIGNMENT_MATCH:
                    alignedQuery_.push_back(query_[qPos++]);
                    alignedTarget_.push_back(target_[tPos++]);
                    break;

                case CigarOperationType::DELETION:
                    alignedQuery_.push_back('-');
                    alignedTarget_.push_back(target_[tPos++]);
                    break;

                case CigarOperationType::INSERTION:
                case CigarOperationType::SOFT_CLIP:
                    alignedQuery_.push_back(query_[qPos++]);
                    alignedTarget_.push_back('-');
                    break;

                case CigarOperationType::HARD_CLIP:
                case CigarOperationType::PADDING:
                case CigarOperationType::REFERENCE_SKIP:
                    throw std::runtime_error(std::string{"unsupported CIGAR op encountered: "} +
                                             op.Char());

                case CigarOperationType::UNKNOWN_OP:
                    throw std::runtime_error("unknown CIGAR op encountered");
            }
        }
    }
}

BandedChainAlignment::BandedChainAlignment(const BandedChainAlignConfig& config, const char* target,
                                           const size_t targetLen, const char* query,
                                           const size_t queryLen, const PacBio::Data::Cigar& cigar)
    : BandedChainAlignment(config, std::string{target, targetLen}, std::string{query, queryLen},
                           cigar)
{
}

float BandedChainAlignment::Identity() const
{
    assert(alignedQuery_.length() == alignedTarget_.length());

    size_t numMatches = 0;
    const size_t len = alignedQuery_.length();
    for (size_t i = 0; i < len; ++i) {
        if (alignedQuery_.at(i) == alignedTarget_.at(i)) ++numMatches;
    }
    return (100.0f * static_cast<float>(numMatches) / len);
}

int64_t BandedChainAlignment::Score() const
{
    using PacBio::Data::CigarOperationType;

    int64_t score = 0;
    const size_t numOps = cigar_.size();
    for (size_t i = 0; i < numOps; ++i) {
        const auto& op = cigar_.at(i);
        switch (op.Type()) {
            case CigarOperationType::SEQUENCE_MATCH:
                score += (config_.matchScore_ * op.Length());
                break;

            case CigarOperationType::SEQUENCE_MISMATCH:
                score += (config_.mismatchPenalty_ * op.Length());
                break;

            case CigarOperationType::INSERTION:
            case CigarOperationType::DELETION:
                // do not penalize end gaps
                if (i != 0 && i != numOps - 1)
                    score +=
                        (config_.gapOpenPenalty_ + config_.gapExtendPenalty_ * (op.Length() - 1));
                break;

            case CigarOperationType::ALIGNMENT_MATCH:
            case CigarOperationType::HARD_CLIP:
            case CigarOperationType::PADDING:
            case CigarOperationType::REFERENCE_SKIP:
            case CigarOperationType::SOFT_CLIP:
                throw std::runtime_error(std::string{"unexpected CIGAR op encountered: "} +
                                         op.Char());

            case CigarOperationType::UNKNOWN_OP:
                throw std::runtime_error("unknown CIGAR op encountered");
        }
    }
    return score;
}

// ------------------------
// BandedChainAlignConfig
// ------------------------

BandedChainAlignConfig BandedChainAlignConfig::Default()
{
    // match, mismatch, gapOpen, gapExtend, band
    return BandedChainAlignConfig{2.0f, -1.0f, -2.0f, -1.0f, 15};
}

// ------------------------
// API 'free functions'
// ------------------------

BandedChainAlignment BandedChainAlign(const char* target, const size_t targetLen, const char* query,
                                      const size_t queryLen,
                                      const std::vector<PacBio::Align::Seed>& seeds,
                                      const BandedChainAlignConfig& config)
{
    Internal::BandedChainAlignerImpl impl{config};
    return impl.Align(target, targetLen, query, queryLen, seeds);
}

}  // namespace Align
}  // namespace PacBio
