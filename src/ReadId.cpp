// Author: Lance Hepler

#include <iostream>
#include <sstream>
#include <string>

#include <pacbio/data/ReadId.h>

using std::ostream;

namespace PacBio {
namespace Data {

ReadId::operator std::string() const
{
    std::ostringstream ss;
    ss << *this;
    return ss.str();
}

ostream& operator<<(ostream& os, const ReadId& id)
{
    os << *(id.MovieName) << '/' << id.HoleNumber;

    if (id.ZmwInterval) {
        os << '/' << id.ZmwInterval->Left() << '_' << id.ZmwInterval->Right();
    }

    return os;
}

}  // namespace Data
}  // namespace PacBio
