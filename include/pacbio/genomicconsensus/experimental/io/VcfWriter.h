// Author: Derek Barnett

#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <pacbio/genomicconsensus/experimental/Settings.h>
#include <pacbio/genomicconsensus/experimental/Variant.h>
#include <pacbio/genomicconsensus/experimental/io/FileProducer.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

class VcfWriter
{
public:
    VcfWriter(const Settings& settings, const std::vector<ReferenceWindow>& refWindows);

    void WriteVariant(const Variant& variant);
    void WriteVariants(const std::vector<Variant>& variants);

private:
    void WriteLine(const std::string& line);

private:
    FileProducer file_;
    std::ofstream out_;
};

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
