// Author: Derek Barnett

#include <iostream>

#include <pacbio/align/PairwiseAlignment.h>
#include <pacbio/genomicconsensus/Input.h>
#include <pacbio/genomicconsensus/Output.h>
#include <pacbio/genomicconsensus/Variant.h>
#include <pacbio/genomicconsensus/arrow/Arrow.h>

using Arrow = PacBio::GenomicConsensus::Arrow;
using Input = PacBio::GenomicConsensus::Input;
using Output = PacBio::GenomicConsensus::Output;
using PairwiseAlignment = PacBio::Align::PairwiseAlignment;
using Settings = PacBio::GenomicConsensus::Settings;
using Variant = PacBio::GenomicConsensus::Variant;

int main(int argc, char* argv[])
{
    if (argc < 3) {
        std::cerr << "usage: arrow <input BAM> <input FASTA>" << std::endl;
        return EXIT_FAILURE;
    }

    Settings settings;
    settings.inputFilename = argv[1];
    settings.referenceFilename = argv[2];

    Input input{settings};
    Output output{settings};
    for (const auto& window : input.ReferenceWindows()) {
        // TODO: some quick cutout for no coverage

        //
        // We call consensus on the enlarged window and then map back
        // to the reference and clip the consensus at the implied
        // bounds.  This seems to be more reliable thank cutting the
        // consensus bluntly.
        //
        const auto enlargedWindow = input.EnlargedWindow(window);
        const auto refSeq = input.ReferenceInWindow(enlargedWindow);
        const auto cssAndVariants =
            Arrow::ConsensusAndVariantsForWindow(input, enlargedWindow, refSeq, settings);

        //
        // Restrict the consensus and variants to the reference window.
        //
        const std::unique_ptr<PairwiseAlignment> ga{
            PacBio::Align::Align(refSeq, cssAndVariants.css.sequence)};

        const auto targetPositions = PacBio::Align::TargetToQueryPositions(*ga);
        const auto cssStart = targetPositions.at(window.Start() - enlargedWindow.Start());
        const auto cssEnd = targetPositions.at(window.End() - enlargedWindow.Start());

        const auto cssSeq = std::string{cssAndVariants.css.sequence.begin() + cssStart,
                                        cssAndVariants.css.sequence.begin() + cssEnd};

        const auto cssQv = std::vector<uint8_t>{cssAndVariants.css.confidence.begin() + cssStart,
                                                cssAndVariants.css.confidence.begin() + cssEnd};

        std::vector<Variant> variants;
        for (const auto& v : cssAndVariants.variants) {
            if (window.Start() <= v.refStart && v.refStart < window.End()) {
                variants.push_back(v);
            }
        }

        output.AddResult({{window, std::move(cssSeq), std::move(cssQv)}, std::move(variants)});
    }
}
