// Author: Lance Hepler

#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <pacbio/data/ChemistryMapping.h>
#include <pacbio/data/ChemistryTriple.h>
#include <pacbio/exception/CCSExceptions.h>
#include <pbcopper/utility/FileUtils.h>

using std::string;

namespace PacBio {
namespace Data {

using namespace PacBio::Utility;
using namespace PacBio::Exception;

ChemistryMapping::ChemistryMapping(const std::string& mappingXML)
{
    using boost::property_tree::ptree;
    using boost::property_tree::xml_parser::read_xml;

    ptree pt;

    if (!FileExists(mappingXML)) throw BadMappingXMLException("File does not exist: " + mappingXML);

    read_xml(mappingXML, pt);

    try {
        mapping[ChemistryTriple::Null()] =
            pt.get<string>("MappingTable.DefaultSequencingChemistry");

        BOOST_FOREACH (ptree::value_type& v, pt.get_child("MappingTable")) {
            if (v.first.compare("Mapping") == 0) {
                ChemistryTriple entry(v.second.get<string>("BindingKit"),
                                      v.second.get<string>("SequencingKit"),
                                      v.second.get<string>("SoftwareVersion"));
                mapping[entry] = v.second.get<string>("SequencingChemistry");
            }
        }
    } catch (...) {
        throw BadMappingXMLException("Could not parse mapping xml!");
    }
}

string ChemistryMapping::MapTriple(const ChemistryTriple& triple, const std::string& fallback) const
{
    try {
        return mapping.at(triple);
    } catch (const std::out_of_range& e) {
        if (fallback.empty()) throw;
        return fallback;
    }
}

}  // namespace Data
}  // namespace PacBio
