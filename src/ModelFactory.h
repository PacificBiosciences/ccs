// Copyright (c) 2011-2016, Pacific Biosciences of California, Inc.
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
