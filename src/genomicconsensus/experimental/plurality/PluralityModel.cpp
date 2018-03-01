// Author: Derek Barnett

#include <pacbio/genomicconsensus/experimental/plurality/PluralityModel.h>

#include <pacbio/genomicconsensus/experimental/Input.h>
#include <pacbio/genomicconsensus/experimental/ReferenceWindow.h>
#include <pacbio/genomicconsensus/experimental/Settings.h>
#include <pacbio/genomicconsensus/experimental/WorkChunk.h>
#include <pacbio/genomicconsensus/experimental/plurality/Plurality.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

WindowResult PluralityModel::ConsensusAndVariantsForWindow(const Input& input,
                                                           const ReferenceWindow& window,
                                                           std::string refSeq,
                                                           const Settings& settings)
{
    return Plurality::ConsensusAndVariantsForWindow(input, window, std::move(refSeq), settings);
}

WindowResult PluralityModel::ProcessChunk(const WorkChunk& chunk, const Settings& settings)
{
    const Input input{settings};
    return ConsensusAndVariantsForWindow(input, chunk.window, input.ReferenceInWindow(chunk.window),
                                         settings);
}

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
