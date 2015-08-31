
#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace PacBio {
namespace Consensus {

struct Read
{
    Read(const std::string& name,
         const std::string& seq,
         const std::vector<uint8_t>& cov,
         const std::string& model);

    std::string Name;
    std::string Seq;
    std::vector<uint8_t> Cov;
    std::string Model;
};

enum struct StrandEnum : uint8_t
{
    FORWARD,
    REVERSE
};

struct MappedRead : public Read
{
    MappedRead(const Read& read,
               StrandEnum strand,
               size_t templateStart,
               size_t templateEnd,
               bool pinStart = false,
               bool pinEnd = false);

    StrandEnum Strand;
    size_t TemplateStart;
    size_t TemplateEnd;
    bool PinStart;
    bool PinEnd;
};

} // namespace Consensus
} // namespace PacBio
