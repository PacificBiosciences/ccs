// Copyright (c) 2015-2016, Pacific Biosciences of California, Inc.
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
#include <set>
#include <string>

#include <boost/optional.hpp>

#include <pacbio/data/Read.h>
#include <pacbio/exception/ModelError.h>

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
    static bool Register(const std::string& name, std::unique_ptr<ModelCreator>&& ctor);
    static boost::optional<std::string> Resolve(const std::string& name);
    static std::set<std::string> SupportedModels();

private:
    static std::map<std::string, std::unique_ptr<ModelCreator>>& CreatorTable();
};

// The concrete form of ModelCreator, which registers a compiled-in
//   model with the ModelFactory and implements the aforementioned
//   Create method for instantiating a concrete model given an SNR
template <typename T>
class ModelCreatorImpl : public ModelCreator
{
public:
    ModelCreatorImpl<T>() {}
    ModelCreatorImpl<T>(const std::set<std::string>& names)
    {
        for (const std::string& name : names)
            if (!ModelFactory::Register(name + "::Compiled",
                                        std::unique_ptr<ModelCreator>(new ModelCreatorImpl<T>())))
                throw PacBio::Exception::DuplicateModel(name);
    }

    virtual std::unique_ptr<ModelConfig> Create(const SNR& snr) const
    {
        return std::unique_ptr<ModelConfig>(new T(snr));
    }
};

// An accessor to a global parameter for overriding the model
boost::optional<std::string>& ModelOverride();

#define REGISTER_MODEL(cls) \
private:                    \
    static const ModelCreatorImpl<cls> creator_

#define REGISTER_MODEL_IMPL(cls) const ModelCreatorImpl<cls> cls::creator_(cls::Names())

}  // namespace Consensus
}  // namespace PacBio
