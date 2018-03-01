// Author: Derek Barnett

#pragma once

#include <string>
#include <vector>

#include <pacbio/genomicconsensus/Consensus.h>
#include <pacbio/genomicconsensus/Variant.h>

namespace PacBio {
namespace GenomicConsensus {

struct WindowResult
{
    Consensus css;
    std::vector<Variant> variants;
};

}  // namespace GenomicConsensus
}  // namespace PacBio