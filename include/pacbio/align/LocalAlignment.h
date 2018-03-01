// Author: Derek Barnett

//
// SIMD local (Smith-Waterman) alignment score
//

#pragma once

#include <string>
#include <vector>

namespace PacBio {
namespace Align {

class LocalAlignment
{
public:
    LocalAlignment(const int32_t targetBegin, const int32_t targetEnd, const int32_t queryBegin,
                   const int32_t queryEnd, const int32_t mismatches, const uint16_t score,
                   std::vector<uint32_t> cigar, std::string cigarString);

    LocalAlignment(const LocalAlignment&) = delete;
    LocalAlignment(LocalAlignment&&) = default;
    LocalAlignment& operator=(const LocalAlignment&) = delete;
    LocalAlignment& operator=(LocalAlignment&&) = default;
    ~LocalAlignment(void) = default;

public:
    int32_t TargetBegin(void) const { return targetBegin_; }
    int32_t TargetEnd(void) const { return targetEnd_; }
    int32_t QueryBegin(void) const { return queryBegin_; }
    int32_t QueryEnd(void) const { return queryEnd_; }
    int32_t NumMismatches(void) const { return mismatches_; }
    uint16_t Score(void) const { return score_; }
    std::vector<uint32_t> Cigar(void) const { return cigar_; }
    std::string CigarString(void) const { return cigarString_; }

private:
    int32_t targetBegin_;
    int32_t targetEnd_;
    int32_t queryBegin_;
    int32_t queryEnd_;
    int32_t mismatches_;
    uint16_t score_;
    std::vector<uint32_t> cigar_;
    std::string cigarString_;
};

struct LocalAlignConfig
{
public:
    uint8_t MatchScore;
    uint8_t MismatchPenalty;
    uint8_t GapOpenPenalty;
    uint8_t GapExtendPenalty;

public:
    static LocalAlignConfig Default(void);
};

///
/// \brief LocalAlign
///
/// \param target
/// \param query
/// \param config
///
/// \return
///
LocalAlignment LocalAlign(const std::string& target, const std::string& query,
                          const LocalAlignConfig& config = LocalAlignConfig::Default());

///
/// \brief LocalAlign
///
/// \param target
/// \param queries
/// \param config
///
/// \return
///
std::vector<LocalAlignment> LocalAlign(
    const std::string& target, const std::vector<std::string>& queries,
    const LocalAlignConfig& config = LocalAlignConfig::Default());

}  // namespace Align
}  // namespace PacBio
