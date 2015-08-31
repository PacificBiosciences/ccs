
#pragma once

#include <array>
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <pacbio/consensus/ParameterTable.h>

namespace PacBio {
namespace Consensus {
namespace detail {

extern uint8_t TranslationTable[256];

} // namespace detail

struct TemplatePosition
{
    char Base;
    double Match;
    double Branch;
    double Stick;
    double Deletion;
};

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
    virtual std::vector<TemplatePosition> Populate(const std::string& tpl) const = 0;
    virtual double BaseEmissionPr(char from, char to) const = 0;
    virtual double CovEmissionPr(MoveType move, uint8_t cov) const = 0;

protected:
    template<typename T>
    static bool Register()
    {
        auto& tbl = ParameterTable::Default_();
        tbl.tbl_[T::Name()] = [] (const SNR& snr)
        {
            return std::unique_ptr<ModelConfig>(new T(snr));
        };
        return true;
    }
};

} // namespace Consensus
} // namespace PacBio
