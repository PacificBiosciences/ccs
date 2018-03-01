// Author: Lance Hepler

#pragma once

#include <map>
#include <stdexcept>
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/optional.hpp>

#include <pacbio/data/IntervalTree.h>

namespace PacBio {
namespace CCS {

class Whitelist
{
public:
    Whitelist(const std::string& spec) : all(false), globalZmws(boost::none)
    {
        // if we're all or *:*, then set all and break out
        if (spec == "all" || spec == "*:*") {
            all = true;
            return;
        }

        std::vector<std::string> mspecs;
        boost::split(mspecs, spec, boost::is_any_of(";"));

        for (const auto& mspec : mspecs) {
            // no craziness policy
            if (mspec == "all" || mspec == "*:*" || globalZmws)
                throw std::invalid_argument("invalid whitelist specification");

            std::vector<std::string> parts;
            boost::split(parts, mspec, boost::is_any_of(":"));

            // only 1 part, it's a ZMW range
            if (parts.size() == 1) {
                if (movieZmws.empty()) {
                    globalZmws = PacBio::Data::IntervalTree::FromString(parts[0]);
                    continue;
                }
            }
            // two parts, but *:_? again, it's just a ZMW range
            else if (parts.size() == 2 && parts[0] == "*") {
                if (movieZmws.empty()) {
                    globalZmws = PacBio::Data::IntervalTree::FromString(parts[1]);
                    continue;
                }
            }
            // two parts, either _:_ or _:*?
            //   either grab a range from the movie or everything, respectively
            else if (parts.size() == 2 && movieZmws.find(parts[0]) == movieZmws.end()) {
                if (parts[1] == "*")
                    movieZmws[parts[0]] = boost::none;
                else
                    movieZmws[parts[0]] = PacBio::Data::IntervalTree::FromString(parts[1]);

                continue;
            }

            // anything else is bad, including resetting any range
            throw std::invalid_argument("invalid whitelist specification");
        }
    }

    bool Contains(const std::string& movieName, int32_t holeNumber) const
    {
        if (all) return true;

        if (globalZmws) return globalZmws->Contains(holeNumber);

        auto it = movieZmws.find(movieName);
        if (it != movieZmws.end()) return !it->second || it->second->Contains(holeNumber);

        return false;
    }

private:
    bool all;
    boost::optional<PacBio::Data::IntervalTree> globalZmws;
    std::map<std::string, boost::optional<PacBio::Data::IntervalTree>> movieZmws;
};

}  // namespace CCS
}  // namespace PacBio
