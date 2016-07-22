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

using ModelError = PacBio::Exception::ModelError;

std::map<std::string, ModelFormCreator*>& ModelFormFactory::CreatorTable()
{
    static std::map<std::string, ModelFormCreator*> tbl;
    return tbl;
}

bool ModelFormFactory::LoadModel(const std::string& path)
{
    using boost::property_tree::ptree;
    using boost::property_tree::read_json;

    ptree pt;

    try {
        read_json(path, pt);
        // verify we're looking at consensus model parameters
        std::string version = pt.get<std::string>("ConsensusModelVersion");
        if (version != "3.0.0") return false;
    } catch (boost::property_tree::ptree_error) {
        return false;
    }

    const std::string chemistry = pt.get<std::string>("ChemistryName");
    const std::string form = pt.get<std::string>("ModelForm");

    const auto tbl = CreatorTable();
    const auto it = tbl.find(form);

    if (it == tbl.end()) return false;

    const std::string name = chemistry + "::" + form + "::" + "FromFile";

    try {
        return ModelFactory::Register(name, it->second->LoadParams(pt));
    } catch (ModelError& e) {
        return false;
    }
}

bool ModelFormFactory::Register(const std::string& form, ModelFormCreator* ctor)
{
    return CreatorTable().insert(std::make_pair(form, ctor)).second;
}
}
}
