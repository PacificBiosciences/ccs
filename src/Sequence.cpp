// Author: Lance Hepler

#include <array>
#include <stdexcept>
#include <string>

using namespace std::literals::string_literals;  // for std::operator ""s

#include <pacbio/data/Sequence.h>

namespace PacBio {
namespace Data {

char Complement(const char base)
{
    constexpr const std::array<char, 256> lookupTable{
        {/*   0 -   7: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /*   8 -  15: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /*  16 -  23: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /*  24 -  31: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /*  32 -  39: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /*  40 -  47: */ 0,   0,   0,   0,   0,   '-', 0,   0,
         /*  48 -  55: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /*  56 -  63: */ 0,   0,   0,   0,   0,   0,   0,   0,

         /*  64 -  71: */ 0,   'T', 'V', 'G', 'H', 0,   0,   'C',
         /*  72 -  79: */ 'D', 0,   0,   'M', 0,   'K', 'N', 0,
         /*  80 -  87: */ 0,   0,   'Y', 'S', 'A', 0,   'B', 'W',
         /*  88 -  95: */ 0,   'R', 0,   0,   0,   0,   0,   0,

         /*  96 - 103: */ 0,   't', 'v', 'g', 'h', 0,   0,   'c',
         /* 104 - 111: */ 'd', 0,   0,   'm', 0,   'k', 'n', 0,
         /* 112 - 119: */ 0,   0,   'y', 's', 'a', 0,   'b', 'w',
         /* 120 - 127: */ 0,   'r', 0,   0,   0,   0,   0,   0,

         /* 128 - 135: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 136 - 143: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 144 - 151: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 152 - 159: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 160 - 167: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 168 - 175: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 176 - 183: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 184 - 191: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 192 - 199: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 200 - 207: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 208 - 215: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 216 - 223: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 224 - 231: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 232 - 239: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 240 - 247: */ 0,   0,   0,   0,   0,   0,   0,   0,
         /* 248 - 255: */ 0,   0,   0,   0,   0,   0,   0,   0}};

    const char result = lookupTable[static_cast<uint8_t>(base)];

    if (result == 0) throw std::invalid_argument(base + " is an invalid base!"s);

    return result;
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

}  // namespace Data
}  // namespace PacBio
