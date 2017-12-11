// Copyright (c) 2016-2017, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: Lance Hepler

#pragma once

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
        return std::unique_ptr<ModelCreator>(new T(pt));
    }
};

#define REGISTER_MODELFORM_IMPL(MODEL)                                   \
    void Init##MODEL()                                                   \
    {                                                                    \
        using namespace MODEL;                                           \
        static const ModelFormCreatorImpl<MODEL##ModelCreator> creator_{ \
            MODEL##ModelCreator::Form()};                                \
    }

}  // namespace Consensus
}  // namespace PacBio
