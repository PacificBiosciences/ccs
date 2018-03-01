// Author: Lance Hepler

#pragma once

#include <iostream>
#include <memory>
#include <string>

#include <boost/optional.hpp>

#include <pacbio/data/Interval.h>

namespace PacBio {
namespace Data {

struct ReadId
{
    std::shared_ptr<std::string> MovieName;
    size_t HoleNumber;
    boost::optional<Interval> ZmwInterval;

    ReadId(const std::shared_ptr<std::string>& movieName, size_t holeNumber)
        : MovieName(movieName), HoleNumber{holeNumber}, ZmwInterval{boost::none}
    {
    }

    ReadId(const std::shared_ptr<std::string>& movieName, size_t holeNumber,
           const Interval& interval)
        : MovieName(movieName), HoleNumber{holeNumber}, ZmwInterval(interval)
    {
    }

    operator std::string() const;
    friend std::ostream& operator<<(std::ostream&, const ReadId&);
};

}  // namespace Data
}  // namespace PacBio
