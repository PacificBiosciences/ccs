// Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.
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

#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <pacbio/ccs/ChemistryMapping.h>
#include <pacbio/ccs/ChemistryTriple.h>
#include <pacbio/ccs/Exceptions.h>
#include <pacbio/ccs/Utility.h>

using namespace std;

namespace PacBio {
namespace CCS {

ChemistryMapping::ChemistryMapping(const std::string& mappingXML)
{
    using boost::property_tree::ptree;
    using boost::property_tree::xml_parser::read_xml;

    ptree pt;

    if (!FileExists(mappingXML))
        throw BadMappingXMLException("File does not exist: " + mappingXML);

    read_xml(mappingXML, pt);

    try
    {
        mapping[ChemistryTriple::Null()] = pt.get<string>("MappingTable.DefaultSequencingChemistry");

        BOOST_FOREACH(ptree::value_type &v, pt.get_child("MappingTable"))
        {
            if (v.first.compare("Mapping") == 0)
            {
                ChemistryTriple entry(v.second.get<string>("BindingKit"),
                                      v.second.get<string>("SequencingKit"),
                                      v.second.get<string>("SoftwareVersion"));
                mapping[entry] = v.second.get<string>("SequencingChemistry");
            }
        }
    }
    catch (...)
    {
        throw BadMappingXMLException("Could not parse mapping xml!");
    }
}

string ChemistryMapping::MapTriple(const ChemistryTriple& triple,
                                   const std::string& fallback) const
{
    try
    {
        return mapping.at(triple);
    }
    catch (const out_of_range& e)
    {
        if (fallback.empty())
            throw;
        return fallback;
    }
}

} // namespace CCS
} // namespace PacBio
