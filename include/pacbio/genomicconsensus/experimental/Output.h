// Author: Derek Barnett

#pragma once

#include <map>
#include <memory>
#include <string>

#include <pacbio/genomicconsensus/experimental/Consensus.h>
#include <pacbio/genomicconsensus/experimental/Settings.h>
#include <pacbio/genomicconsensus/experimental/Variant.h>
#include <pacbio/genomicconsensus/experimental/WindowResult.h>
#include <pacbio/genomicconsensus/experimental/io/FastaWriter.h>
#include <pacbio/genomicconsensus/experimental/io/FastqWriter.h>
#include <pacbio/genomicconsensus/experimental/io/GffWriter.h>
#include <pacbio/genomicconsensus/experimental/io/VcfWriter.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

class Output
{
public:
    explicit Output(const Settings& settings);

public:
    void AddResult(WindowResult result);

private:
    void MaybeFlushContig(const std::string& refName);

private:
    Settings settings_;

    // writers
    std::unique_ptr<FastaWriter> fasta_;
    std::unique_ptr<FastqWriter> fastq_;
    std::unique_ptr<GffWriter> gff_;
    std::unique_ptr<VcfWriter> vcf_;

    // per-reference data
    std::map<std::string, ReferenceWindow> refWindows_;
    std::map<std::string, uint32_t> expectedBasesPerRef_;
    std::map<std::string, uint32_t> processedBasesPerRef_;
    std::map<std::string, std::vector<Consensus>> consensiPerRef_;
    std::map<std::string, std::vector<Variant>> variantsPerRef_;
};

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
