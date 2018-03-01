// Author: Derek Barnett

#include <pacbio/genomicconsensus/experimental/poa/PoaModel.h>

#include <pacbio/genomicconsensus/experimental/ReferenceWindow.h>
#include <pacbio/genomicconsensus/experimental/Settings.h>
#include <pacbio/genomicconsensus/experimental/poa/Poa.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

WindowResult PoaModel::ConsensusAndVariantsFromWindow(
    const Input& /*input*/, const std::vector<PacBio::BAM::BamRecord>& reads,
    const ReferenceWindow& window, const std::string& refSeq, const Settings& settings) const
{
    return Poa::ConsensusAndVariantsForAlignments(window, refSeq, reads, settings);
}

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
