// Author: Derek Barnett

#include <pacbio/genomicconsensus/experimental/io/FastaWriter.h>

#include <stdexcept>
#include <vector>

using namespace std::literals::string_literals;

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

namespace {

std::vector<std::string> WrapToColumns(const std::string& seq, const size_t width = 70)
{
    std::vector<std::string> result;
    size_t pos = 0;
    while (pos < seq.size()) {
        result.push_back(seq.substr(pos, width));
        pos += width;
    }
    return result;
}

}  // namespace anonymous

FastaWriter::FastaWriter(const Settings& settings)
    : file_{settings.fastaFilename}, out_{file_.tempFilename_}
{
    //  targetFilename, targetFilename + ".tmp"
    if (!out_) {
        throw std::runtime_error("could not open "s + file_.targetFilename_ + " for writing"s);
    }
}

void FastaWriter::Write(const std::string& header, const std::string& sequence)
{
    WriteLine('>' + header);
    for (const auto& line : WrapToColumns(sequence))
        WriteLine(line);
}

void FastaWriter::WriteLine(const std::string& line) { out_ << line << '\n'; }

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
