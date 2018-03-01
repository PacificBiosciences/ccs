// Author: Lance Hepler

#include <iostream>
#include <sstream>
#include <string>

#include <pacbio/data/Interval.h>

namespace PacBio {
namespace Data {

Interval::operator std::string() const
{
    std::ostringstream ss;
    ss << *this;
    return ss.str();
}

std::ostream& operator<<(std::ostream& os, const Interval& interval)
{
    return os << '[' << interval.Left() << ", " << interval.Right() << ')';
}

}  // namespace Data
}  // namespace PacBio
