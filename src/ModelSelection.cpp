// Author: Lance Hepler

#include <mutex>
#include <set>
#include <string>
#include <vector>

#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <boost/optional.hpp>

#include <pacbio/consensus/ModelSelection.h>

#include "ModelFactory.h"
#include "ModelFormFactory.h"
#include "ModelNaming.h"

namespace PacBio {
namespace Consensus {

std::set<std::string> SupportedModels() { return ModelFactory::SupportedModels(); }

std::set<std::string> SupportedChemistries()
{
    std::set<std::string> chems;
    for (const auto& model : SupportedModels())
        chems.insert(model.substr(0, model.find_first_of(':')));
    return chems;
}

bool OverrideModel(const std::string& name)
{
    boost::optional<std::string> model(boost::none);

    if (!(model = ModelFactory::Resolve(name))) return false;

    ModelOverride() = *model;
    return true;
}

bool UnOverrideModel()
{
    ModelOverride() = boost::none;
    return true;
}

bool LoadModelFromFile(const std::string& path, const ModelOrigin origin)
{
    struct stat st;
    if (stat(path.c_str(), &st) != 0 || !S_ISREG(st.st_mode)) return false;
    return ModelFormFactory::LoadModel(path, origin);
}

boost::optional<size_t> LoadModelsFromDirectory(const std::string& dirPath,
                                                const ModelOrigin origin, const bool strict)
{
    static std::mutex m;

    struct stat st;
    if (stat(dirPath.c_str(), &st) != 0) return boost::none;
    if (!S_ISDIR(st.st_mode)) return boost::none;

    DIR* dp = opendir(dirPath.c_str());
    if (dp == nullptr) return boost::none;

    size_t nModels = 0, dot;
    {  // Lock down this block to prevent multiple calls to readdir()
        std::lock_guard<std::mutex> lock(m);

        // iterate through .json files in directory,
        //   loading any into ModelFactory
        bool ret = true;
        struct dirent* ep;
        while ((ep = readdir(dp)) != nullptr) {
            std::string path = dirPath + '/' + ep->d_name;
            if ((dot = path.find_last_of('.')) == std::string::npos || path.substr(dot) != ".json")
                continue;
            if (!(ret &= (stat(path.c_str(), &st) == 0))) break;
            if (S_ISREG(st.st_mode)) {
                if ((ret &= LoadModelFromFile(path, origin))) {
                    ++nModels;
                } else if (strict) {
                    closedir(dp);
                    return boost::none;
                } else
                    break;
            }
        }

        closedir(dp);
    }  // End of lock_guard block

    return boost::make_optional(nModels);
}

size_t LoadModels(const std::string& path)
{
    static const std::string origin = "FromFile";
    struct stat st;
    if (stat(path.c_str(), &st) != 0)
        return 0;
    else if (S_ISDIR(st.st_mode))
        return LoadModelsFromDirectory(path, origin, false).get_value_or(0);
    else if (S_ISREG(st.st_mode))
        return LoadModelFromFile(path, origin) ? 1 : 0;
    return 0;
}
}
}
