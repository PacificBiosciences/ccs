// Author: Derek Barnett

#pragma once

#include <string>

#include <pacbio/genomicconsensus/experimental/IConsensusModel.h>
#include <pacbio/genomicconsensus/experimental/WindowResult.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

class Input;
struct ReferenceWindow;
struct Settings;
struct WorkChunk;

///
/// \brief The PluralityModel class
///
class PluralityModel : public IConsensusModel
{
public:
    PluralityModel() = default;
    ~PluralityModel() = default;

    ///
    /// \brief ConsensusAndVariantsForWindow
    /// \param input
    /// \param window
    /// \param refSeq
    /// \param settings
    /// \return
    ///
    WindowResult ConsensusAndVariantsForWindow(const Input& input, const ReferenceWindow& window,
                                               std::string refSeq, const Settings& settings);

    ///
    /// \brief ProcessChunk
    /// \param chunk
    /// \param settings
    /// \return
    ///
    WindowResult ProcessChunk(const WorkChunk& chunk, const Settings& settings);
};

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
