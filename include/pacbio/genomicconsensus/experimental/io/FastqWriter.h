// Author: Derek Barnett

#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <pacbio/genomicconsensus/experimental/Settings.h>
#include <pacbio/genomicconsensus/experimental/io/FileProducer.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

class FastqWriter
{
public:
    FastqWriter(const Settings& settings);
    void Write(const std::string& header, const std::string& sequence,
               const std::vector<uint8_t>& qualities);

private:
    void WriteLine(const std::string& line);

private:
    FileProducer file_;
    std::ofstream out_;
};

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
