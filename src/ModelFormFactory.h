// Author: Lance Hepler

#pragma once

#include "UnanimityInternalConfig.h"

#include <map>
#include <memory>
#include <string>

#include <boost/property_tree/ptree.hpp>

#include "ModelFactory.h"
#include "ModelNaming.h"

namespace PacBio {
namespace Consensus {

// this pattern is based on
// http://blog.fourthwoods.com/2011/06/04/factory-design-pattern-in-c/

// An abstract class defining a single abstract method, LoadParams, for
//   parameterizing a model form yielding a ModelCreator that can be
//   used to instantiate a concrete model given an SNR (see ModelCreator)
class ModelFormCreator
{
public:
    virtual ~ModelFormCreator() {}
    virtual std::unique_ptr<ModelCreator> LoadParams(
        const boost::property_tree::ptree& pt) const = 0;
};

// A static factory class that holds onto all available model forms,
//   with methods to: get the static map of ModelFormCreators accessible
//   by their form name; to register a form within said map; and to load
//   a model parameter file, find its model form, and insert a
//   parameterized model into the ModelFactory where it is discoverable
//   by the rest of the library
class ModelFormFactory
{
public:
    static bool LoadModel(const std::string& path, const ModelOrigin origin);
    static bool Register(ModelForm form, ModelFormCreator* ctor);

private:
    static std::map<ModelForm, ModelFormCreator*>& CreatorTable();
};

// The concrete version of the ModelFormCreator, which registers
//   a model form with the ModelFormFactory, and implements the
//   aforementioned abstract method LoadParams to yield a ModelCreator
//   given a property_tree full of parameters
template <typename T>
class ModelFormCreatorImpl : public ModelFormCreator
{
public:
    ModelFormCreatorImpl<T>(const ModelForm form)
    {
        if (!ModelFormFactory::Register(form, this))
            throw std::runtime_error("duplicate model form inserted into form factory!");
    }

    virtual std::unique_ptr<ModelCreator> LoadParams(const boost::property_tree::ptree& pt) const
    {
        return std::make_unique<T>(pt);
    }
};

#define REGISTER_MODELFORM_IMPL(MODEL)                                   \
    UNANIMITY_PRIVATE_API void Init##MODEL()                             \
    {                                                                    \
        using namespace MODEL;                                           \
        static const ModelFormCreatorImpl<MODEL##ModelCreator> creator_{ \
            MODEL##ModelCreator::Form()};                                \
    }

}  // namespace Consensus
}  // namespace PacBio
