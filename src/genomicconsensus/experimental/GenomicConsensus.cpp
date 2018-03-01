// Author: Derek Barnett

#include <pacbio/genomicconsensus/experimental/GenomicConsensus.h>

#include <pbbam/IndexedFastaReader.h>

#include <pacbio/genomicconsensus/experimental/Consensus.h>
#include <pacbio/genomicconsensus/experimental/ConsensusModelFactory.h>
#include <pacbio/genomicconsensus/experimental/Settings.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

WindowResult Process(const WorkChunk& chunk, const Settings& settings)
{
    const auto& window = chunk.window;
    if (!chunk.hasCoverage) {
        // quick skip
        PacBio::BAM::IndexedFastaReader fasta{settings.referenceFilename};
        const auto refSeq = fasta.Subsequence(window.name, window.Start(), window.End());
        return WindowResult{Consensus::NoCallConsensus(settings.noCallStyle, window, refSeq),
                            std::vector<Variant>{}};
    } else {
        auto model = ConsensusModelFactory::Create(settings.mode);
        return model->ProcessChunk(chunk, settings);
    }
}

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
