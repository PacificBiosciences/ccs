// Author: Armin Töpfer

#pragma once

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include <pacbio/data/StrandType.h>

namespace PacBio {
namespace Data {

/// Stores nucleotide-wise signal to noise ratios.
struct SNR
{
    double A;
    double C;
    double G;
    double T;

    SNR(double a, double c, double g, double t);
    SNR(const std::vector<double>& snrs);
    SNR(const std::vector<float>& snrs);

    SNR(const double (&snrs)[4]) : A{snrs[0]}, C{snrs[1]}, G{snrs[2]}, T{snrs[3]} {}

    operator std::vector<float>() const
    {
        std::vector<float> snr = {static_cast<float>(A), static_cast<float>(C),
                                  static_cast<float>(G), static_cast<float>(T)};
        return snr;
    }

    inline const double& operator[](const size_t i) const
    {
        if (i == 0) return A;
        if (i == 1) return C;
        if (i == 2) return G;
        if (i == 3) return T;
        throw std::invalid_argument("SNR out of bounds!");
    }

    inline double& operator[](const size_t i)
    {
        // casting away const when underlying object is non-const, is well-defined
        return const_cast<double&>(static_cast<const SNR&>(*this)[i]);
    }

    inline bool operator==(const SNR& other) const
    {
        return other.A == A && other.C == C && other.G == G && other.T == T;
    }

    inline bool operator!=(const SNR& other) const { return !(*this == other); }

    inline double Minimum(void) const { return std::min(std::min(A, C), std::min(G, T)); }
};

SNR ClampSNR(const SNR& val, const SNR& min, const SNR& max);

/// A Read contains the name, sequence, covariates, SNR, and associated model.
struct Read
{
    Read(const std::string& name, const std::string& seq, const std::vector<uint8_t>& ipd,
         const std::vector<uint8_t>& pw, const SNR& snr, std::string model);
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

/// A MappedRead extends Read by the strand information and template anchoring
/// positions.
struct MappedRead : public Read
{
    MappedRead(const Read& read, StrandType strand, size_t templateStart, size_t templateEnd,
               bool pinStart = false, bool pinEnd = false);
    MappedRead(const MappedRead& read) = default;
    MappedRead(MappedRead&& read) = default;

    StrandType Strand;
    size_t TemplateStart;
    size_t TemplateEnd;
    bool PinStart;
    bool PinEnd;
};

std::ostream& operator<<(std::ostream&, const MappedRead&);

}  // namespace Data
}  // namespace PacBio
