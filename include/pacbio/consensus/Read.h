
#pragma once

#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace PacBio {
namespace Consensus {

struct Read
{
    Read(const std::string& name, const std::string& seq, const std::vector<uint8_t>& cov,
         const std::string& model);
    Read(const std::string& name, const std::string& seq, const std::string& model);
    Read(const Read& read) = default;
    Read(Read&& read) = default;

    std::string Name;
    std::string Seq;
    std::vector<uint8_t> Cov;
    std::string Model;

    inline size_t Length() const { return Seq.length(); }
};

enum struct StrandEnum : uint8_t
{
    FORWARD,
    REVERSE,
    UNMAPPED
};

struct MappedRead : public Read
{
    MappedRead(const Read& read, StrandEnum strand, size_t templateStart, size_t templateEnd,
               bool pinStart = false, bool pinEnd = false);
    MappedRead(const MappedRead& read) = default;
    MappedRead(MappedRead&& read) = default;

    StrandEnum Strand;
    size_t TemplateStart;
    size_t TemplateEnd;
    bool PinStart;
    bool PinEnd;
};

std::ostream& operator<<(std::ostream&, const MappedRead&);

}  // namespace Consensus
}  // namespace PacBio
