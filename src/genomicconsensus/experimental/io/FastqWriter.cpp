// Author: Derek Barnett

#include <pacbio/genomicconsensus/experimental/io/FastqWriter.h>

#include <stdexcept>

#include <pbbam/QualityValues.h>

using namespace std::literals::string_literals;

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

FastqWriter::FastqWriter(const Settings& settings)
    : file_{settings.fastqFilename}, out_{file_.tempFilename_}
{
    //  targetFilename, targetFilename + ".tmp"
    if (!out_) {
        throw std::runtime_error("could not open "s + file_.targetFilename_ + " for writing"s);
    }
}

void FastqWriter::Write(const std::string& header, const std::string& sequence,
                        const std::vector<uint8_t>& qualities)
{
    WriteLine('@' + header);
    WriteLine(sequence);
    WriteLine("+");
    WriteLine(PacBio::BAM::QualityValues(qualities).Fastq());
}

void FastqWriter::WriteLine(const std::string& line) { out_ << line << '\n'; }

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
