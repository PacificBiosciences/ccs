// Author: Lance Hepler

#pragma once

#include <map>
#include <string>

#include <pacbio/data/ChemistryTriple.h>

namespace PacBio {
namespace Data {

class ChemistryMapping
{
public:
    ChemistryMapping(const std::string& mappingXML);

    std::string MapTriple(const ChemistryTriple& triple, const std::string& fallback = "") const;

private:
    struct ChemistryComparer
    {
        bool operator()(const ChemistryTriple& a, const ChemistryTriple& b) const
        {
            if (a.BindingKit < b.BindingKit or a.SequencingKit < b.SequencingKit or
                a.MajorVersion < b.MajorVersion or a.MinorVersion < b.MinorVersion) {
                return true;
            }
            return false;
        }
    };

    std::map<ChemistryTriple, std::string, ChemistryComparer> mapping;
};

}  // namespace Data
}  // namespace PacBio
