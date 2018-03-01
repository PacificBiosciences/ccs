// Author: Derek Barnett

#include <pacbio/genomicconsensus/experimental/arrow/ArrowModel.h>

#include <pbbam/BamRecord.h>

#include <pacbio/genomicconsensus/experimental/arrow/Arrow.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

WindowResult ArrowModel::ConsensusAndVariantsFromWindow(
    const Input& /*input*/, const std::vector<PacBio::BAM::BamRecord>& reads,
    const ReferenceWindow& window, const std::string& refSeq, const Settings& settings) const
{
    Consensus css;
    boost::optional<std::vector<uint8_t>> effectiveSiteCoverage;
    if (settings.reportEffectiveCoverage) {
        std::vector<PacBio::BAM::BamRecord> readsUsed;
        css = Arrow::ConsensusForAlignments(window, refSeq, reads, settings, &readsUsed);
        effectiveSiteCoverage = Arrow::CoverageInWindow(window, readsUsed);
    } else
        css = Arrow::ConsensusForAlignments(window, refSeq, reads, settings, nullptr);

    const auto siteCoverage = Arrow::CoverageInWindow(window, reads);
    return Arrow::VariantsFromConsensus(window, refSeq, css, siteCoverage, effectiveSiteCoverage,
                                        settings);
}

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
