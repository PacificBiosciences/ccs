// Author: Lance Hepler

#include <cassert>
#include <cstdint>

#include <pacbio/consensus/ModelConfig.h>

namespace PacBio {
namespace Consensus {

std::ostream& operator<<(std::ostream& out, const TemplatePosition& pos)
{
    return out << "TemplatePosition(" << pos.Base << ", " << pos.Match << ", " << pos.Branch << ", "
               << pos.Stick << ", " << pos.Deletion << ')';
}

}  // namespace Consensus
}  // namespace PacBio
