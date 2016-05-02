
#pragma once

#include <map>
#include <memory>
#include <set>
#include <string>

namespace PacBio {
namespace Consensus {

// forward declarations
struct SNR;
class ModelConfig;

// this pattern is based on
// http://blog.fourthwoods.com/2011/06/04/factory-design-pattern-in-c/
class ModelCreator
{
public:
    ModelCreator(const std::set<std::string>& names);
    virtual ~ModelCreator() {}
    virtual std::unique_ptr<ModelConfig> Create(const SNR&) const = 0;
};

template <typename T>
class ModelCreatorImpl : public ModelCreator
{
public:
    ModelCreatorImpl<T>(const std::set<std::string>& names) : ModelCreator(names) {}
    virtual std::unique_ptr<ModelConfig> Create(const SNR& snr) const
    {
        return std::unique_ptr<ModelConfig>(new T(snr));
    }
};

class ModelFactory
{
public:
    static std::unique_ptr<ModelConfig> Create(const std::string& name, const SNR&);
    static bool Register(const std::string& name, ModelCreator* ctor);
    static std::set<std::string> SupportedChemistries();

private:
    static std::map<std::string, ModelCreator*>& CreatorTable();
};

#define REGISTER_MODEL(cls) \
private:                    \
    static const ModelCreatorImpl<cls> creator_

#define REGISTER_MODEL_IMPL(cls) const ModelCreatorImpl<cls> cls::creator_(cls::Names())

}  // namespace Consensus
}  // namespace PacBio
