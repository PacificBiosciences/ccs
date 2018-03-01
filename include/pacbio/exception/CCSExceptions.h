// Author: Lance Hepler

#pragma once

#include <cstddef>
#include <stdexcept>

namespace PacBio {
namespace Exception {

class BadMappingXMLException : public std::runtime_error
{
public:
    BadMappingXMLException(const std::string& msg) : std::runtime_error(msg) {}
};

class BadChemistryTriple : public std::runtime_error
{
public:
    BadChemistryTriple(const std::string& msg) : std::runtime_error(msg) {}
};

}  // namespace Exception
}  // namespace PacBio
