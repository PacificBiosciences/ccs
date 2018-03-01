// Author: Lance Hepler

#pragma once

#include <stdexcept>
#include <string>

namespace PacBio {
namespace Exception {

class ModelError : public std::runtime_error
{
public:
    ModelError(const std::string& msg) : std::runtime_error(msg) {}
};

class ChemistryNotFound : public ModelError
{
public:
    ChemistryNotFound(const std::string& name)
        : ModelError(std::string("chemistry not found: '") + name + "'")
    {
    }
};

class DuplicateModel : public ModelError
{
public:
    DuplicateModel(const std::string& name) : ModelError("duplicate model: '" + name + "'") {}
};

class MalformedModelFile : public ModelError
{
public:
    MalformedModelFile() : ModelError("malformed model!") {}
};

}  // namespace Exception
}  // namespace PacBio
