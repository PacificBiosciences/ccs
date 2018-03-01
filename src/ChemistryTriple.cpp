// Author: Lance Hepler

#include <boost/lexical_cast.hpp>
#include <boost/xpressive/xpressive.hpp>

#include <pacbio/data/ChemistryTriple.h>
#include <pacbio/exception/CCSExceptions.h>

namespace PacBio {
namespace Data {

using namespace PacBio::Exception;

ChemistryTriple::ChemistryTriple(const std::string& bindingKit, const std::string& sequencingKit,
                                 const std::string& changeListID)
{
    if (!SetValues(bindingKit, sequencingKit, changeListID)) {
        throw BadChemistryTriple("unparsable ChemistryTriple(" + bindingKit + ", " + sequencingKit +
                                 ", " + changeListID + ")");
    }
}

bool ChemistryTriple::SetValues(const std::string& bindingKit, const std::string& sequencingKit,
                                const std::string& changeListID)
{
    using namespace boost::xpressive;

    using boost::lexical_cast;

    try {
        BindingKit = lexical_cast<unsigned>(bindingKit);
        SequencingKit = lexical_cast<unsigned>(sequencingKit);

        smatch what;
        // sregex::compile("^(\\d+)\\.(\\d+)");
        sregex re = bos >> (s1 = +_d) >> '.' >> (s2 = +_d);

        if (regex_search(changeListID.begin(), changeListID.end(), what, re)) {
            MajorVersion = lexical_cast<unsigned>(what[1]);
            MinorVersion = lexical_cast<unsigned>(what[2]);
            return true;
        }
    } catch (const boost::bad_lexical_cast& e) {
        return false;
    }

    return false;
}

}  // namespace Data
}  // namespace PacBio
