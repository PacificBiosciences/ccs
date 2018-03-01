// Author: Derek Barnett

#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <tuple>

#include <boost/optional.hpp>

namespace PacBio {
namespace GenomicConsensus {

struct Variant
{
public:
    std::string refName;
    size_t refStart;
    size_t refEnd;
    std::string refSeq;
    std::string readSeq1;
    char refPrev;
    char readPrev;

    boost::optional<std::string> readSeq2 = boost::none;
    boost::optional<size_t> frequency1 = boost::none;
    boost::optional<size_t> frequency2 = boost::none;
    boost::optional<size_t> coverage = boost::none;
    boost::optional<size_t> confidence = boost::none;
    boost::optional<std::map<std::string, std::string> > annotations = boost::none;

public:
    Variant(const std::string& _refName, const size_t _refStart, const size_t _refEnd,
            const std::string& _refSeq, const std::string& _readSeq1, const char _refPrev,
            const char _readPrev)
        : refName{_refName}
        , refStart{_refStart}
        , refEnd{_refEnd}
        , refSeq{_refSeq}
        , readSeq1{_readSeq1}
        , refPrev{_refPrev}
        , readPrev{_readPrev}
    {
    }

    Variant() = default;
    Variant(const Variant&) = default;
    Variant(Variant&&) = default;
    Variant& operator=(const Variant&) = default;
    Variant& operator=(Variant&&) = default;
    ~Variant() = default;

public:
    void Annotate(const std::string& key, const std::string& value)
    {
        annotations->insert(std::make_pair(key, value));
    }

    bool IsHeterozygous() const { return (bool)readSeq2; }

    bool operator<(const Variant& other) const
    {
        return std::tie(refName, refStart, refEnd, readSeq1) <
               std::tie(other.refName, other.refStart, other.refEnd, other.readSeq1);
    }
};

}  // namespace GenomicConsensus
}  // namespace PacBio
