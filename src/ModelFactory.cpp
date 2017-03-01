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

#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/optional.hpp>

#include <pacbio/consensus/ModelConfig.h>
#include <pacbio/exception/StateError.h>

#include "ModelFactory.h"

namespace PacBio {
namespace Consensus {

// fwd decl from ModelSelection.cpp
boost::optional<size_t> LoadModelsFromDirectory(const std::string& path, const ModelOrigin origin,
                                                bool strict);

namespace {

using ChemistryNotFound = PacBio::Exception::ChemistryNotFound;

size_t Count(const std::string& str, const std::string& delim)
{
    size_t count = 0;
    size_t start = 0;

    while ((start = str.find(delim, start)) != std::string::npos) {
        ++count;
        start += delim.length();
    }

    return count;
}

void LoadBundleModels()
{
    static bool updatesLoaded = false;
    if (!updatesLoaded) {
        if (const char* pth = getenv("PB_CHEMISTRY_BUNDLE_DIR")) {
            if (!LoadModelsFromDirectory(std::string(pth) + "/arrow", ModelOrigin::BUNDLED, true))
                throw Exception::ModelError(
                    std::string("unable to load arrow model updates from: ") + pth);
            updatesLoaded = true;
        }
    }
}
}

std::unique_ptr<ModelConfig> ModelFactory::Create(const std::string& name, const SNR& snr)
{
    // Load update bundle models before we create anything
    LoadBundleModels();

    boost::optional<std::string> model(boost::none);

    if (!(model = ModelOverride()))
        if (!(model = Resolve(name))) throw ChemistryNotFound(name);

    const auto& tbl = CreatorTable();
    const auto it = tbl.find(*model);

    if (it == tbl.end()) throw ChemistryNotFound(name);

    return it->second->Create(snr);
}

std::unique_ptr<ModelConfig> ModelFactory::Create(const PacBio::Data::Read& read)
{
    return Create(read.Model, read.SignalToNoise);
}

bool ModelFactory::Register(const ModelName& name, std::unique_ptr<ModelCreator>&& ctor)
{
    return CreatorTable()
        .insert(std::make_pair(name, std::forward<std::unique_ptr<ModelCreator>>(ctor)))
        .second;
}

boost::optional<std::string> ModelFactory::Resolve(const std::string& name)
{
    const std::vector<std::string> forms = ModelForm::Preferences();
    const std::vector<std::string> origins = ModelOrigin::Preferences();
    const auto& tbl = CreatorTable();
    const size_t nParts = Count(name, "::") + 1;

    if (nParts == 3) {
        if (tbl.find(name) != tbl.end()) return name;
    }

    else if (nParts == 2) {
        for (const auto& origin : origins) {
            const std::string model = name + "::" + origin;
            if (tbl.find(model) != tbl.end()) return model;
        }
    }

    else if (nParts == 1) {
        for (const auto& form : forms) {
            for (const auto& origin : origins) {
                const std::string model = name + "::" + form + "::" + origin;
                if (tbl.find(model) != tbl.end()) return model;
            }
        }
    }

    return boost::none;
}

std::set<std::string> ModelFactory::SupportedModels()
{
    // Load update bundle models before we report anything
    LoadBundleModels();

    const auto& tbl = CreatorTable();
    std::set<std::string> result;
    for (const auto& item : tbl)
        result.insert(item.first);
    return result;
}

std::map<ModelName, std::unique_ptr<ModelCreator>>& ModelFactory::CreatorTable()
{
    static std::map<ModelName, std::unique_ptr<ModelCreator>> tbl;
    return tbl;
}

boost::optional<std::string>& ModelOverride()
{
    static boost::optional<std::string> ovr = boost::none;
    return ovr;
}

}  // Consensus
}  // PacBio
