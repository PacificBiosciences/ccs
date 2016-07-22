// Copyright (c) 2016, Pacific Biosciences of California, Inc.
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
#include <string>
#include <vector>

#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <pacbio/consensus/ModelSelection.h>

#include "ModelFactory.h"
#include "ModelFormFactory.h"

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

bool LoadModelFromFile(const std::string& path)
{
    struct stat st;
    if (stat(path.c_str(), &st) != 0 || !S_ISREG(st.st_mode)) return false;
    return ModelFormFactory::LoadModel(path);
}

size_t LoadModelsFromDirectory(const std::string& dirPath)
{
    struct stat st;
    if (stat(dirPath.c_str(), &st) != 0) return false;
    if (!S_ISDIR(st.st_mode)) return false;

    DIR* dp = opendir(dirPath.c_str());
    if (dp == nullptr) return false;

    // iterate through .json files in directory,
    //   loading any into ModelFactory
    bool ret = true;
    size_t nModels = 0;
    struct dirent* ep;
    while ((ep = readdir(dp)) != nullptr) {
        std::string path = dirPath + '/' + ep->d_name;
        if (path.substr(path.find_last_of('.')) != ".json") continue;
        if (!(ret &= (stat(path.c_str(), &st) == 0))) break;
        if (S_ISREG(st.st_mode)) {
            if ((ret &= LoadModelFromFile(path)))
                ++nModels;
            else
                break;
        }
    }

    closedir(dp);

    return nModels;
}

size_t LoadModels(const std::string& path)
{
    struct stat st;
    if (stat(path.c_str(), &st) != 0)
        return 0;
    else if (S_ISDIR(st.st_mode))
        return LoadModelsFromDirectory(path);
    else if (S_ISREG(st.st_mode))
        return LoadModelFromFile(path) ? 1 : 0;
    return 0;
}
}
}
