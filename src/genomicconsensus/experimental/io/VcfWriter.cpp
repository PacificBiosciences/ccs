// Author: Derek Barnett

#include <pacbio/genomicconsensus/experimental/io/VcfWriter.h>

#include <chrono>
#include <sstream>
#include <stdexcept>

#include <pbbam/DataSet.h>

#include <pacbio/UnanimityVersion.h>

using namespace std::literals::string_literals;

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

VcfWriter::VcfWriter(const Settings& settings, const std::vector<ReferenceWindow>& refWindows)
    : file_{settings.vcfFilename}, out_{file_.tempFilename_}
{
    //  targetFilename, targetFilename + ".tmp"
    if (!out_) {
        throw std::runtime_error("could not open "s + file_.targetFilename_ + " for writing"s);
    }

    WriteLine("##fileformat=VCFv4.3");
    WriteLine("##fileDate=" + PacBio::BAM::ToIso8601(std::chrono::system_clock::now()));
    WriteLine("##source=" + PacBio::UnanimityVersion());
    WriteLine("##reference=file://" + settings.referenceFilename);

    for (const auto& ref : refWindows)
        WriteLine("##contig=<ID=" + ref.name + ",length=" + std::to_string(ref.Length()) + ">");

    WriteLine("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
}

void VcfWriter::WriteLine(const std::string& line) { out_ << line << '\n'; }

void VcfWriter::WriteVariant(const Variant& v)
{
    // "VCF lite"
    auto pos = v.refStart;
    std::string ref;
    std::string alt;

    // indel
    if (v.refSeq.empty() || v.readSeq1.empty() || (v.IsHeterozygous() && v.readSeq2->empty())) {
        ref = v.refPrev + v.refSeq;
        alt = v.readPrev + v.readSeq1;
        if (v.IsHeterozygous()) {
            alt.push_back(',');
            alt += v.readPrev;
            alt += v.readSeq2.get();
        }
    }

    // substitution
    else {
        ++pos;
        ref = v.refSeq;
        if (v.IsHeterozygous()) {
            if (v.refSeq == v.readSeq1)
                alt = v.readSeq2.get();
            else if (v.refSeq == v.readSeq2.get())
                alt = v.readSeq1;
            else {
                alt = v.readSeq1;
                alt.push_back(',');
                alt += v.readSeq2.get();
            }
        } else
            alt = v.readSeq1;
    }

    std::stringstream line;
    line << v.refName << '\t' << pos << "\t.\t" << ref << '\t' << alt << '\t' << v.confidence.get()
         << "\tPASS";
    WriteLine(line.str());
}

void VcfWriter::WriteVariants(const std::vector<Variant>& variants)
{
    for (const auto& v : variants)
        WriteVariant(v);
}

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
