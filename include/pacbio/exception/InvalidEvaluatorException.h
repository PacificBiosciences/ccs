// Author: Armin Töpfer

#pragma once

#include <stdexcept>
#include <string>

namespace PacBio {
namespace Exception {

class InvalidEvaluatorException : public std::runtime_error
{
public:
    InvalidEvaluatorException(const std::string& msg) : std::runtime_error(msg) {}
};

}  // namespace Exception
}  // namespace PacBio