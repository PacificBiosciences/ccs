// Author: Lance Hepler

#include <map>
#include <string>
#include <utility>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <pacbio/exception/ModelError.h>

#include "ModelFactory.h"
#include "ModelFormFactory.h"

namespace PacBio {
namespace Consensus {

std::map<ModelForm, ModelFormCreator*>& ModelFormFactory::CreatorTable()
{
    static std::map<ModelForm, ModelFormCreator*> tbl;
    return tbl;
}

bool ModelFormFactory::LoadModel(const std::string& path, const ModelOrigin origin)
{
    using boost::property_tree::ptree;
    using boost::property_tree::read_json;

    ptree pt;

    try {
        read_json(path, pt);
        // verify we're looking at consensus model parameters
        std::string version = pt.get<std::string>("ConsensusModelVersion");
        if (version != "3.0.0") return false;

        const std::string chemistry = pt.get<std::string>("ChemistryName");
        const ModelForm form(pt.get<std::string>("ModelForm"));

        const auto tbl = CreatorTable();
        const auto it = tbl.find(form);

        if (it == tbl.end()) return false;

        const ModelName name(chemistry, form, origin);

        return ModelFactory::Register(name, it->second->LoadParams(pt));
    } catch (boost::property_tree::ptree_error&) {
    } catch (Exception::ModelNamingError&) {
    } catch (Exception::ModelError& e) {
    }
    return false;
}

bool ModelFormFactory::Register(const ModelForm form, ModelFormCreator* ctor)
{
    return CreatorTable().insert(std::make_pair(form, ctor)).second;
}
}
}
