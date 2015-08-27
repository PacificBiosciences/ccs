
#pragma once

#include <array>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace PacBio {
namespace Consensus {

typedef std::array<double, 4> SNR;

struct TemplatePosition
{
    char Base;
    double Match;
    double Branch;
    double Stick;
    double Deletion;
};

enum MoveType
{
    MATCH    = 0,
    BRANCH   = 1,
    STICK    = 2,
    DELETION = 3 // never used for covariate
};

class ModelConfig
{
public:
    virtual std::vector<TemplatePosition> Populate(
            const std::string& tpl) const = 0;
    virtual double CovariatePr(MoveType move, uint8_t cov) const = 0;
    virtual double MiscallPr(char from, char to) const = 0;
};

class ParameterTable
{
public:
    bool Contains(const std::string& name) const;
    std::unique_ptr<ModelConfig> At(const std::string& name, const SNR& snr) const;

    static const ParameterTable& Default();

private:
    static ParameterTable& Default_();

    std::map<std::string, std::function<ModelConfig(const SNR&)>> tbl_;

    friend class ModelConfig;
};

} // namespace Consensus
} // namespace PacBio
