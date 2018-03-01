// Author: Derek Barnett

#include <pacbio/genomicconsensus/experimental/Output.h>

#include <iostream>

#include <pbcopper/logging/Logging.h>

#include <pacbio/genomicconsensus/experimental/Input.h>
#include <pacbio/genomicconsensus/experimental/Workflow.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

Output::Output(const Settings& settings) : settings_{settings}
{
    const Input input{settings_};
    const auto refWindows = input.ReferenceWindows(false);

    for (const auto& window : refWindows) {
        refWindows_[window.name] = window;
        processedBasesPerRef_[window.name] = 0;
        expectedBasesPerRef_[window.name] = window.Length();
        consensiPerRef_[window.name] = {};
        variantsPerRef_[window.name] = {};
    }

    // init writers from settings
    if (!settings_.fastaFilename.empty()) fasta_ = std::make_unique<FastaWriter>(settings_);
    if (!settings_.fastqFilename.empty()) fastq_ = std::make_unique<FastqWriter>(settings_);
    if (!settings_.gffFilename.empty()) gff_ = std::make_unique<GffWriter>(settings_, refWindows);
    if (!settings_.vcfFilename.empty()) vcf_ = std::make_unique<VcfWriter>(settings_, refWindows);
}

void Output::AddResult(WindowResult result)
{
    PBLOG_INFO << "Adding result for " << result.css.window;

    const auto window = result.css.window;
    consensiPerRef_[window.name].emplace_back(std::move(result.css));
    variantsPerRef_[window.name];  // append list
    processedBasesPerRef_[window.name] += window.Length();
    MaybeFlushContig(window.name);
}

void Output::MaybeFlushContig(const std::string& refName)
{
    const auto basesProcessed = processedBasesPerRef_[refName];
    const auto expectedBases = expectedBasesPerRef_[refName];
    if (basesProcessed == expectedBases) {
        // This contig is done, so we can dump to file & delete data
        if (gff_ || vcf_) {
            // sort & write variants
            auto& variants = variantsPerRef_[refName];
            std::sort(variants.begin(), variants.end());
            if (gff_) gff_->WriteVariants(variants);
            if (vcf_) vcf_->WriteVariants(variants);
        }

        // TODO: this doesn't actually free up memory
        variantsPerRef_[refName].clear();

        //
        // If the user asked to analyze a window or a set of
        // windows, we output a FAST[AQ] contig per analyzed
        // window.  Otherwise we output a fasta contig per
        // reference contig.
        //
        // We try to be intelligent about naming the output
        // contigs, to include window information where applicable.
        //
        for (const auto& window : Workflow::EnumerateWindows(refName, settings_)) {
            std::string recordName = window.name;
            if (window != refWindows_[window.name]) {
                recordName +=
                    '_' + std::to_string(window.Start()) + '_' + std::to_string(window.End());
            }

            std::string algoName;
            switch (settings_.mode) {
                case ConsensusMode::ARROW:
                    algoName = "arrow";
                    break;
                case ConsensusMode::PLURALITY:
                    algoName = "plurality";
                    break;
                case ConsensusMode::POA:
                    algoName = "poa";
                    break;
                default:
                    throw std::runtime_error("unknown consensus mode");
            }
            std::string cssName = recordName + '|' + algoName;

            std::vector<Consensus> consensiInThisWindow;
            for (const auto& consensus : consensiPerRef_[refName]) {
                if (Overlap(consensus.window, window)) consensiInThisWindow.push_back(consensus);
            }
            auto css = Consensus::Join(consensiInThisWindow);

            if (fasta_) fasta_->Write(cssName, css.sequence);
            if (fastq_) fastq_->Write(cssName, css.sequence, css.confidence);
        }

        // TODO: this doesn't actually free up memory
        consensiPerRef_[refName].clear();
    }
}

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
