// Author: Derek Barnett

#pragma once

#include <string>
#include <vector>

#include <pbbam/BamRecord.h>

#include <pacbio/genomicconsensus/experimental/IPoaModel.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

struct ReferenceWindow;
struct Settings;

///
/// \brief The PoaModel class
///
class PoaModel : public IPoaModel
{
public:
    PoaModel() = default;
    ~PoaModel() = default;

    ///
    /// \brief ConsensusAndVariantsFromWindow
    /// \param input
    /// \param reads
    /// \param window
    /// \param refSeq
    /// \param settings
    /// \return
    ///
    WindowResult ConsensusAndVariantsFromWindow(const Input& /*input*/,
                                                const std::vector<PacBio::BAM::BamRecord>& reads,
                                                const ReferenceWindow& window,
                                                const std::string& refSeq,
                                                const Settings& settings) const;
};

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
