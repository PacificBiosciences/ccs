// Author: Armin Töpfer

#pragma once

#include <string>

namespace PacBio {
namespace Data {

char Complement(char base);
std::string Complement(const std::string& input);
std::string Reverse(const std::string& input);
std::string ReverseComplement(const std::string& input);

}  // namespace Data
}  // namespace PacBio
