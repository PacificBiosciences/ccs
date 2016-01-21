
#pragma once

#include <array>
#include <cstdint>
#include <map>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace PacBio {
namespace Consensus {
namespace detail {

extern uint8_t TranslationTable[256];

}  // namespace detail

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
};

struct TemplatePosition
{
    char Base;
    double Match;
    double Branch;
    double Stick;
    double Deletion;
};

std::ostream& operator<<(std::ostream&, const TemplatePosition&);

enum struct MoveType : uint8_t
{
    MATCH = 0,
    BRANCH = 1,
    STICK = 2,
    DELETION = 3  // never used for covariate
};

class ModelConfig
{
public:
    virtual ~ModelConfig() {}
    virtual std::vector<TemplatePosition> Populate(const std::string& tpl) const = 0;
    virtual double BaseEmissionPr(MoveType move, char from, char to) const = 0;
    virtual double CovEmissionPr(MoveType move, uint8_t cov, const char from, const char to) const = 0;
    // folded into CovEmissionPr for now:
    //   virtual double CounterWeight() const = 0;
    virtual double UndoCounterWeights(size_t nEmissions) const = 0;
};

}  // namespace Consensus
}  // namespace PacBio
