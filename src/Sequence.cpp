
#include <stdexcept>
#include <string>

#include <pacbio/consensus/Sequence.h>

namespace PacBio {
namespace Consensus {

char Complement(char base)
{
    switch (base) {
        case 'A':
            return 'T';
        case 'a':
            return 't';
        case 'C':
            return 'G';
        case 'c':
            return 'g';
        case 'G':
            return 'C';
        case 'g':
            return 'c';
        case 'T':
            return 'A';
        case 't':
            return 'a';
        case '-':
            return '-';
        default:
            throw std::invalid_argument("invalid base!");
    }
}

std::string Complement(const std::string& input)
{
    std::string output;
    output.reserve(input.length());
    for (const char b : input)
        output.push_back(Complement(b));
    return output;
}

std::string Reverse(const std::string& input)
{
    std::string output;
    output.reserve(input.length());
    for (auto it = input.crbegin(); it != input.crend(); ++it)
        output.push_back(*it);
    return output;
}

std::string ReverseComplement(const std::string& input)
{
    std::string output;
    output.reserve(input.length());
    for (auto it = input.crbegin(); it != input.crend(); ++it)
        output.push_back(Complement(*it));
    return output;
}

}  // namespace Consensus
}  // namespace PacBio
