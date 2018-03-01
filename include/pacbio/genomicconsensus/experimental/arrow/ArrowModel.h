// Author: Derek Barnett

#pragma once

#include <pacbio/genomicconsensus/experimental/IPoaModel.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

///
/// \brief The ArrowModel class
///
class ArrowModel : public IPoaModel
{
public:
    ArrowModel() = default;
    ~ArrowModel() = default;

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
