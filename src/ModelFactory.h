// Author: Lance Hepler

#pragma once

#include "UnanimityInternalConfig.h"

#include <map>
#include <memory>
#include <set>
#include <string>

#include <boost/optional.hpp>

#include <pacbio/data/Read.h>
#include <pacbio/exception/ModelError.h>

#include "ModelNaming.h"

namespace PacBio {

// forward declarations
namespace Data {
struct SNR;
}

namespace Consensus {

// forward declarations
class ModelConfig;

using SNR = PacBio::Data::SNR;

// this pattern is based on
// http://blog.fourthwoods.com/2011/06/04/factory-design-pattern-in-c/

// An abstract class with a single abstract method, Create, for
//   instantiating a concrete model given the discriminative SNR
class ModelCreator
{
public:
    virtual ~ModelCreator() {}
    virtual std::unique_ptr<ModelConfig> Create(const SNR&) const = 0;
};

// A static class containing the map of parameterized models that need
//   only SNR to become concrete, with methods to create such models,
//   register models, resolve models, and list available models
class ModelFactory
{
public:
    static std::unique_ptr<ModelConfig> Create(const std::string& name, const SNR&);
    static std::unique_ptr<ModelConfig> Create(const PacBio::Data::Read& read);
    static bool Register(const ModelName& name, std::unique_ptr<ModelCreator>&& ctor);
    static boost::optional<std::string> Resolve(const std::string& name);
    static std::set<std::string> SupportedModels();

private:
    static std::map<ModelName, std::unique_ptr<ModelCreator>>& CreatorTable();
};

// The concrete form of ModelCreator, which registers a compiled-in
//   model with the ModelFactory and implements the aforementioned
//   Create method for instantiating a concrete model given an SNR
template <typename T>
class ModelCreatorImpl : public ModelCreator
{
public:
    ModelCreatorImpl<T>() {}
    ModelCreatorImpl<T>(const std::set<std::string>& chemistries, const ModelForm form)
    {
        for (const std::string& chemistry : chemistries)
            if (!ModelFactory::Register({chemistry, form, ModelOrigin::COMPILED},
                                        std::make_unique<ModelCreatorImpl<T>>()))
                throw PacBio::Exception::DuplicateModel(chemistry);
    }

    virtual std::unique_ptr<ModelConfig> Create(const SNR& snr) const
    {
        return std::make_unique<T>(snr);
    }
};

// An accessor to a global parameter for overriding the model
boost::optional<std::string>& ModelOverride();

#define REGISTER_MODEL_IMPL(MODEL)                                                          \
    UNANIMITY_PRIVATE_API void Init##MODEL()                                                \
    {                                                                                       \
        using namespace MODEL;                                                              \
        static const ModelCreatorImpl<MODEL##_Model> creator_{MODEL##_Model::Chemistries(), \
                                                              MODEL##_Model::Form()};       \
    }

}  // namespace Consensus
}  // namespace PacBio
