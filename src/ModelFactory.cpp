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
        const char* pth = getenv("SMRT_CHEMISTRY_BUNDLE_DIR");
        if (pth != nullptr && pth[0] != '\0') {
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
    return CreatorTable().emplace(name, std::move(ctor)).second;
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
