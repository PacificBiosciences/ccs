
#pragma once

#include <stdexcept>
#include <string>

namespace PacBio {
namespace Consensus {

class AlphaBetaMismatch : public std::runtime_error
{
public:
    AlphaBetaMismatch() : std::runtime_error("alpha/beta mismatch error!") {}
};

class ChemistryNotFound : public std::runtime_error
{
public:
    ChemistryNotFound(const std::string& name)
        : std::runtime_error(std::string("chemistry not found: '") + name + "'")
    {
    }
};

}  // namespace Consensus
}  // namespace PacBio
