// Author: Derek Barnett

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include <pacbio/genomicconsensus/experimental/NoCallStyle.h>
#include <pacbio/genomicconsensus/experimental/ReferenceWindow.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

///
/// \brief The Consensus struct
///
struct Consensus
{
    ReferenceWindow window;
    std::string sequence;
    std::vector<uint8_t> confidence;

    ///
    /// \brief NoCallConsensus
    /// \param style
    /// \param window
    /// \param refSeq
    /// \return
    ///
    static Consensus NoCallConsensus(const NoCallStyle style, const ReferenceWindow& window,
                                     const std::string& refSeq);

    ///
    /// \brief Join
    /// \param subconsensi
    /// \return
    ///
    static Consensus Join(std::vector<Consensus> subconsensi);
};

inline bool operator<(const Consensus& lhs, const Consensus& rhs)
{
    return lhs.window < rhs.window;
}
inline bool operator==(const Consensus& lhs, const Consensus& rhs)
{
    return lhs.window == rhs.window;
}
inline bool operator!=(const Consensus& lhs, const Consensus& rhs)
{
    return lhs.window != rhs.window;
}

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
