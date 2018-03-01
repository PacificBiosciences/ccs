// Author: Lance Hepler

#include "Recursor.h"

#include <utility>

namespace PacBio {
namespace Consensus {

AbstractRecursor::AbstractRecursor(PacBio::Data::MappedRead mr, const double scoreDiff)
    : read_{std::move(mr)}, scoreDiff_{exp(scoreDiff)}
{
}

}  // namespace Consensus
}  // namespace PacBio
