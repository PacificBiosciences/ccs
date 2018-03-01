// Author: Derek Barnett

#pragma once

#include <string>

namespace PacBio {
namespace GenomicConsensus {

enum class SortingStrategy
{
    LONGEST_AND_STRAND_BALANCED,
    LONGEST,
    SPANNING,
    FILE_ORDER
};

}  // namespace GenomicConsensus
}  // namespace PacBio