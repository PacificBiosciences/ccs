
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

// fwd decl
class AbstractRecursor;
class AbstractTemplate;
struct MappedRead;
struct SNR;

struct TemplatePosition
{
    char Base;
    uint8_t Idx;
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
    virtual std::unique_ptr<AbstractRecursor> CreateRecursor(
        std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr, double scoreDiff) const = 0;
    virtual std::vector<TemplatePosition> Populate(const std::string& tpl) const = 0;
    virtual double SubstitutionRate(uint8_t prev, uint8_t curr) const = 0;
};

}  // namespace Consensus
}  // namespace PacBio
