
#pragma once

#include <array>
#include <cstdint>
#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

namespace PacBio {
namespace Consensus {
namespace detail {

extern uint8_t TranslationTable[256];

} // namespace detail

typedef std::array<double, 4> SNR;

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
    MATCH    = 0,
    BRANCH   = 1,
    STICK    = 2,
    DELETION = 3 // never used for covariate
};

class ModelConfig
{
public:
    virtual ~ModelConfig() { }
    virtual std::vector<TemplatePosition> Populate(const std::string& tpl) const = 0;
    virtual double BaseEmissionPr(char from, char to) const = 0;
    virtual double CovEmissionPr(MoveType move, uint8_t cov) const = 0;
};

} // namespace Consensus
} // namespace PacBio
