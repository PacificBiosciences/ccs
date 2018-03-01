// Author: Lance Hepler

#pragma once

#include <string>

namespace PacBio {
namespace Data {

class ChemistryTriple
{
public:
    unsigned BindingKit;
    unsigned SequencingKit;
    unsigned MajorVersion;
    unsigned MinorVersion;

    ChemistryTriple() : BindingKit{0}, SequencingKit{0}, MajorVersion{0}, MinorVersion{0} {}

    ChemistryTriple(const std::string& bindingKit, const std::string& sequencingKit,
                    const std::string& changeListID);

    static ChemistryTriple Null() { return ChemistryTriple(); }

    bool IsNull() const
    {
        return BindingKit == 0 and SequencingKit == 0 and MajorVersion == 0 and MinorVersion == 0;
    }

    void SetNull()
    {
        BindingKit = 0;
        SequencingKit = 0;
        MajorVersion = 0;
        MinorVersion = 0;
    }

    bool SetValues(const std::string& bindingKit, const std::string& sequencingKit,
                   const std::string& changeListID);
};

}  // namespace Data
}  // namespace PacBio
