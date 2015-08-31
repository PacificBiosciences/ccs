
#pragma once

#include <map>
#include <memory>
#include <string>

namespace PacBio {
namespace Consensus {

typedef std::array<double, 4> SNR;

// forward declaration
class ModelConfig;

class ParameterTable
{
public:
    bool Contains(const std::string& name) const;
    std::unique_ptr<ModelConfig> At(const std::string& name, const SNR& snr) const;

    static const ParameterTable& Default();

private:
    // for Default_()/model registration
    friend class ModelConfig;

    static ParameterTable& Default_();

    std::map<std::string, std::function<std::unique_ptr<ModelConfig>(const SNR&)>> tbl_;
};

} // namespace Consensus
} // namespace PacBio
