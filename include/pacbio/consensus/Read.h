
#pragma once

#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace PacBio {
namespace Consensus {

struct SNR
{
    double A;
    double C;
    double G;
    double T;

    SNR(double a, double c, double g, double t);
    SNR(const std::vector<double>& snrs);

    inline double operator[](const size_t i) const
    {
        if (i == 0) return A;
        if (i == 1) return C;
        if (i == 2) return G;
        if (i == 3) return T;
        throw std::invalid_argument("SNR out of bounds!");
    }

    inline bool operator==(const SNR& other) const
    {
        return other.A == A && other.C == C && other.G == G && other.T == T;
    }

    inline bool operator!=(const SNR& other) const { return !(*this == other); }
};

struct Read
{
    Read(const std::string& name, const std::string& seq, const std::vector<uint8_t>& ipd,
         const std::vector<uint8_t>& pw, const SNR& snr, const std::string& model);
    Read(const Read& read) = default;
    Read(Read&& read) = default;

    std::string Name;
    std::string Seq;
    std::vector<uint8_t> IPD;
    std::vector<uint8_t> PulseWidth;
    SNR SignalToNoise;
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
