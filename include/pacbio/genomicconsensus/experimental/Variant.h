// Author: Derek Barnett

#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <tuple>

#include <boost/optional.hpp>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

///
/// \brief The Variant struct
///
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
    std::map<std::string, std::string> annotations;

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

public:
    ///
    /// \brief Annotate
    /// \param key
    /// \param value
    ///
    void Annotate(const std::string& key, const std::string& value)
    {
        annotations.insert(std::make_pair(key, value));
    }

    ///
    /// \brief IsHeterozygous
    /// \return
    ///
    bool IsHeterozygous() const { return readSeq2.is_initialized() && readSeq1 != readSeq2.get(); }

    bool IsHomozygous() const { return !IsHeterozygous(); }
};

inline bool operator<(const Variant& lhs, const Variant& rhs)
{
    return std::tie(lhs.refName, lhs.refStart, lhs.refEnd, lhs.readSeq1) <
           std::tie(rhs.refName, rhs.refStart, rhs.refEnd, rhs.readSeq1);
}

inline std::string VariantType(const Variant& v)
{
    const auto refLen = v.refSeq.size();
    const auto allele1Len = v.readSeq1.size();
    const auto allele2Len = (v.readSeq2 ? v.readSeq2.get().size() : std::string::npos);

    // ref length is zero
    auto isInsertion = [&]() { return refLen == 0; };

    // either allele length is zero (allele 2 length will be string::npos if missing)
    auto isDeletion = [&]() { return (allele1Len == 0 || allele2Len == 0); };

    // allele1 length is same as reference, and allele2 length is the same as well, if present
    auto isSubstitution = [&]() {
        return ((allele1Len == refLen) &&
                ((allele2Len == std::string::npos) || (allele2Len == refLen)));
    };

    if (isInsertion())
        return "insertion";
    else if (isDeletion())
        return "deletion";
    else if (isSubstitution())
        return "substitution";
    else
        return "variant";
}

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
