
#include <string>

namespace PacBio {
namespace Consensus {

char Complement(char base);
std::string Complement(const std::string& input);
std::string Reverse(const std::string& input);
std::string ReverseComplement(const std::string& input);

}  // namespace Consensus
}  // namespace PacBio
