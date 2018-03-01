// Author: Derek Barnett

#include <pacbio/genomicconsensus/experimental/IPoaModel.h>

#include <pacbio/align/PairwiseAlignment.h>
#include <pacbio/genomicconsensus/experimental/Consensus.h>
#include <pacbio/genomicconsensus/experimental/Filters.h>
#include <pacbio/genomicconsensus/experimental/Input.h>
#include <pacbio/genomicconsensus/experimental/Intervals.h>
#include <pacbio/genomicconsensus/experimental/WorkChunk.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

IPoaModel::~IPoaModel() {}

void IPoaModel::AnnotateVariants(std::vector<Variant>* const variants,
                                 const std::vector<PacBio::BAM::BamRecord>& reads) const
{
    for (auto& v : *variants) {
        std::string annotation;
        for (const auto& read : reads) {
            if (!annotation.empty()) {
                annotation.push_back(',');
                annotation.push_back(' ');
            }
            annotation.append(read.FullName());
        }
        v.Annotate("rows", annotation);
    }
}

void IPoaModel::ClipReadsToWindow(std::vector<BAM::BamRecord>* const reads,
                                  const ReferenceWindow& window) const
{
    const auto winStart = window.Start();
    const auto winEnd = window.End();
    for (auto& read : *reads)
        read.Clip(PacBio::BAM::ClipType::CLIP_TO_REFERENCE, winStart, winEnd);
}

ReferenceWindow IPoaModel::EnlargedWindow(const ReferenceWindow& window, const size_t seqLength,
                                          const size_t overhang) const
{
    const PacBio::Data::Interval refInterval{0, seqLength};
    const auto wStart = window.Start();
    const auto left = ((wStart < overhang) ? 0 : wStart - overhang);
    const auto right = window.End() + overhang;
    return ReferenceWindow{window.name, refInterval.Intersect({left, right})};
}

WindowResult IPoaModel::ProcessChunk(const WorkChunk& chunk, const Settings& settings)
{
    // input reference window
    const Input input{settings};
    const auto& referenceWindow = chunk.window;
    const auto& refName = referenceWindow.name;
    const auto& refSeqLength = input.SequenceLength(refName);

    // enlarged window
    const auto eWindow = EnlargedWindow(referenceWindow, refSeqLength, settings.windowOverhang);
    const auto eStart = eWindow.Start();
    const auto refContig = input.ReferenceInWindow(ReferenceWindow{refName, {0, refSeqLength}});
    const auto refSeqInEnlargedWindow = refContig.substr(eStart, eWindow.Length());

    // CSS/variant calls on enlarged window
    const auto windowResult = ResultForWindow(eWindow, refContig, settings);
    ;

    // restrict CSS/variants to in put window
    auto windowConsensus =
        RestrictedConsensus(windowResult.css, refSeqInEnlargedWindow, referenceWindow);
    auto windowVariants = RestrictedVariants(windowResult.variants, referenceWindow);
    return WindowResult{std::move(windowConsensus), std::move(windowVariants)};
}

WindowResult IPoaModel::ResultForWindow(const ReferenceWindow& refWindow, const std::string& refSeq,
                                        const Settings& settings) const
{
    const auto& winId = refWindow.name;
    const auto& winStart = refWindow.Start();
    const auto& winEnd = refWindow.End();
    const Input input{settings};

    std::vector<Consensus> subconsensi;
    std::vector<Variant> variants;

    // determine intervals for window
    std::vector<PacBio::Data::Interval> allIntervals;
    if (settings.usingFancyChunking) {
        const PacBio::BAM::BamFile bam{settings.inputFilename};
        const PacBio::BAM::PbiRawData index{bam.PacBioIndexFilename()};
        allIntervals = FancyIntervals(index, refWindow, settings);
    } else {
        allIntervals.emplace_back(winStart, winEnd);
    }

    for (const auto& interval : allIntervals) {
        const auto intStart = interval.Left();
        const auto intEnd = interval.Right();
        const auto intRefSeq = refSeq.substr(intStart, interval.Length());
        const auto subWindow = ReferenceWindow{winId, interval};

        auto reads = input.ReadsInWindow(subWindow);
        ClipReadsToWindow(&reads, subWindow);
        FilterAlignments(&reads, settings);

        // determine if this intervals is a "k-spanning" or a "hole"
        const size_t numSpanning = std::count_if(
            reads.begin(), reads.end(), [&subWindow](const PacBio::BAM::BamRecord& read) {
                const auto readStart = static_cast<size_t>(read.ReferenceStart());
                const auto readEnd = static_cast<size_t>(read.ReferenceEnd());
                return readStart <= subWindow.interval.Left() &&
                       subWindow.interval.Right() <= readEnd;
            });

        WindowResult intervalWindow;
        if (numSpanning >= settings.minPoaCoverage) {

            // model-specific CSS/variant calls
            intervalWindow =
                ConsensusAndVariantsFromWindow(input, reads, subWindow, intRefSeq, settings);

            // add variants that pass filters
            if (!intervalWindow.variants.empty()) {
                FilterVariants(&intervalWindow.variants, settings);
                if (settings.annotateGFF) AnnotateVariants(&intervalWindow.variants, reads);
                variants.insert(variants.end(), intervalWindow.variants.begin(),
                                intervalWindow.variants.end());
            }

            // TODO (DB): MaybeDumpEvidence()
        }

        else {
            intervalWindow.css =
                Consensus::NoCallConsensus(NoCallStyle::LOWERCASE_REFERENCE, subWindow, intRefSeq);
        }

        // save interval CSS
        subconsensi.push_back(intervalWindow.css);
    }

    return WindowResult{Consensus::Join(subconsensi), std::move(variants)};
}

Consensus IPoaModel::RestrictedConsensus(const Consensus& enlargedCss, const std::string& refSeq,
                                         const ReferenceWindow& originalWindow) const
{
    const std::unique_ptr<PacBio::Align::PairwiseAlignment> ga{
        PacBio::Align::Align(refSeq, enlargedCss.sequence)};

    const auto targetPositions = PacBio::Align::TargetToQueryPositions(*ga);
    const auto cssStart = targetPositions.at(originalWindow.Start() - enlargedCss.window.Start());
    const auto cssEnd = targetPositions.at(originalWindow.End() - enlargedCss.window.Start());

    const auto cssSeq =
        std::string{enlargedCss.sequence.begin() + cssStart, enlargedCss.sequence.begin() + cssEnd};
    const auto cssQv = std::vector<uint8_t>{enlargedCss.confidence.begin() + cssStart,
                                            enlargedCss.confidence.begin() + cssEnd};

    return Consensus{originalWindow, std::move(cssSeq), std::move(cssQv)};
}

std::vector<Variant> IPoaModel::RestrictedVariants(const std::vector<Variant>& enlargedVariants,
                                                   const ReferenceWindow& originalWindow) const
{
    std::vector<Variant> windowVariants;
    for (const auto& v : enlargedVariants) {
        if (originalWindow.Start() <= v.refStart && v.refStart < originalWindow.End()) {
            windowVariants.push_back(v);
        }
    }
    return windowVariants;
}

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
