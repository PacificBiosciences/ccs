
#pragma once

namespace PacBio {
namespace Consensus {

class AlphaBetaMismatch : public std::runtime_error
{
public:
    AlphaBetaMismatch() : std::runtime_error("alpha/beta mismatch error!")
    { }
};

} // namespace Consensus
} // namespace PacBio
