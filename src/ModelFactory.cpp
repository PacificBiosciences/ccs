
#include <stdexcept>
#include <utility>

#include <pacbio/consensus/Exceptions.h>
#include <pacbio/consensus/ModelConfig.h>

#include "ModelFactory.h"

namespace PacBio {
namespace Consensus {

ModelCreator::ModelCreator(const std::set<std::string>& names)
{
    for (const std::string& name : names)
        if (!ModelFactory::Register(name, this))
            throw std::runtime_error("duplicate model inserted into factory!");
}

std::unique_ptr<ModelConfig> ModelFactory::Create(const std::string& name, const SNR& snr)
{
    const auto it = CreatorTable().find(name);

    if (it == CreatorTable().end()) throw ChemistryNotFound(name);

    return it->second->Create(snr);
}

bool ModelFactory::Register(const std::string& name, ModelCreator* const ctor)
{
    return CreatorTable().insert(std::make_pair(name, ctor)).second;
}

std::set<std::string> ModelFactory::SupportedChemistries()
{
    std::set<std::string> result;
    for (const auto& item : ModelFactory::CreatorTable())
        result.insert(item.first);
    return result;
}

std::map<std::string, ModelCreator*>& ModelFactory::CreatorTable()
{
    static std::map<std::string, ModelCreator*> tbl;
    return tbl;
}

}  // Consensus
}  // PacBio
