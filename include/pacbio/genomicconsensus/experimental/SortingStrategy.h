// Author: Derek Barnett

#pragma once

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

///
/// \brief The SortingStrategy enum
///
enum class SortingStrategy
{
    LONGEST_AND_STRAND_BALANCED,
    LONGEST,
    SPANNING,
    FILE_ORDER
};

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
