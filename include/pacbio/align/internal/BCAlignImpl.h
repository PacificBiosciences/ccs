// Author: Derek Barnett

//
// BandedChainAlignment implementation
//

#pragma once

#include <cassert>

#include <pacbio/align/BandedChainAlignment.h>
#include <pacbio/align/internal/BCAlignBlocks.h>

namespace PacBio {
namespace Align {
namespace Internal {

class BandedChainAlignerImpl
{
public:
    BandedChainAlignerImpl(const BandedChainAlignConfig& config);

public:
    BandedChainAlignment Align(const char* target, const size_t targetLen, const char* query,
                               const size_t queryLen,
                               const std::vector<PacBio::Align::Seed>& seeds);

    void StitchCigars(PacBio::Data::Cigar* global, PacBio::Data::Cigar&& local);

private:
    struct Sequences
    {
        const char* target;
        size_t targetLen;
        const char* query;
        size_t queryLen;
    };

    void AlignGapBlock(const PacBio::Align::Seed& nextSeed);
    void AlignGapBlock(const size_t hLength, const size_t vLength);
    void AlignLastGapBlock(void);

    void AlignSeedBlock(const PacBio::Align::Seed& seed);

    void Initialize(const char* target, const size_t targetLen, const char* query,
                    const size_t queryLen);

    std::vector<PacBio::Align::Seed> MergeSeeds(const std::vector<PacBio::Align::Seed>& seeds);

    BandedChainAlignment Result(void);

private:
    const BandedChainAlignConfig& config_;

    StandardGlobalAlignBlock gapBlock_;
    BandedGlobalAlignBlock seedBlock_;
    PacBio::Data::Cigar globalCigar_;
    int64_t globalScore_;
    size_t gapBlockBeginH_;
    size_t gapBlockBeginV_;
    Sequences sequences_;
};

}  // namespace Internal
}  // namespace Align
}  // namespace PacBio
