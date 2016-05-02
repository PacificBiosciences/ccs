
#include "Recursor.h"

namespace PacBio {
namespace Consensus {

AbstractRecursor::AbstractRecursor(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr,
                                   const double scoreDiff)
    : tpl_{std::forward<std::unique_ptr<AbstractTemplate>>(tpl)}
    , read_{mr}
    , scoreDiff_{exp(scoreDiff)}
{
}

}  // namespace Consensus
}  // namespace PacBio
