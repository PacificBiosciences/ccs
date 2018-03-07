// Author: Derek Barnett

#include <pacbio/genomicconsensus/experimental/io/GffWriter.h>

#include <chrono>
#include <sstream>
#include <stdexcept>

#include <pbbam/DataSet.h>
#include <boost/optional.hpp>

#include <pacbio/UnanimityVersion.h>
#include <pacbio/data/StrandType.h>

using namespace std::literals::string_literals;

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

namespace {

struct GffRecord
{
    // required
    std::string seqId_;
    size_t start_;
    size_t end_;
    std::string type_;

    // optional
    std::string score_;
    std::string strand_;
    std::string phase_;
    std::string source_;
    std::map<std::string, std::string> attributes_;

    GffRecord(std::string seqId, size_t start, size_t end, std::string type,
              std::string score = ".", std::string strand = ".", std::string phase = ".",
              std::string source = ".", std::map<std::string, std::string> attributes = {})
        : seqId_{std::move(seqId)}
        , start_{start}
        , end_{end}
        , type_{std::move(type)}
        , score_{std::move(score)}
        , strand_{std::move(strand)}
        , phase_{std::move(phase)}
        , source_{std::move(source)}
        , attributes_{std::move(attributes)}
    {
    }
};

static GffRecord FromVariant(const Variant& v)
{
    const auto start = v.refStart + (v.refSeq.empty() ? 0 : 1);
    const auto end = (v.refSeq.empty() ? v.refStart : v.refEnd);

    GffRecord gffRecord{v.refName, start, end, VariantType(v)};

    gffRecord.attributes_["reference"] = (v.refSeq.empty() ? "." : v.refSeq);

    gffRecord.attributes_["variantSeq"] = [&]() {
        const auto readSeq1 = (v.readSeq1.empty() ? "." : v.readSeq1);
        if (!v.IsHeterozygous()) return readSeq1;

        // otherwise, format heterozygous
        std::string readSeq2(1, '.');
        if (v.readSeq2 && !v.readSeq2->empty()) readSeq2 = v.readSeq2.get();
        return readSeq1 + '/' + readSeq2;
    }();

    if (v.frequency1) {
        gffRecord.attributes_["frequency"] = [&]() {
            auto freq = std::to_string(v.frequency1.get());
            if (v.IsHeterozygous()) freq += '/' + std::to_string(v.frequency2.get());
            return freq;
        }();
    }

    if (v.coverage) gffRecord.attributes_["coverage"] = std::to_string(v.coverage.get());

    if (v.confidence) gffRecord.attributes_["confidence"] = std::to_string(v.confidence.get());

    // add'l annotations
    for (const auto& att : v.annotations)
        gffRecord.attributes_[att.first] = att.second;

    return gffRecord;
}

}  // anonymous

GffWriter::GffWriter(const Settings& settings, const std::vector<ReferenceWindow>& refWindows)
    : file_{settings.gffFilename}, out_{file_.tempFilename_}
{
    //  targetFilename, targetFilename + ".tmp"
    if (!out_) {
        throw std::runtime_error("could not open "s + file_.targetFilename_ + " for writing"s);
    }

    WriteLine("##gff-version 3");
    WriteLine("##pacbio-variant-version 2.1");
    WriteLine("##date " + PacBio::BAM::ToIso8601(std::chrono::system_clock::now()));
    WriteLine(
        "##feature-ontology"
        " http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.12");
    WriteLine("##source GenomicConsensus " + PacBio::UnanimityVersion());
    WriteLine("##source-commandline " + settings.commandLine);
    WriteLine("##source-alignment-file " + settings.inputFilename);
    WriteLine("##source-reference-file " + settings.referenceFilename);

    for (const auto& ref : refWindows)
        WriteLine("##sequence-region " + ref.name + " 1 " + std::to_string(ref.Length()));
}

void GffWriter::WriteLine(const std::string& line) { out_ << line << '\n'; }

void GffWriter::WriteVariant(const Variant& variant)
{
    const auto gff = FromVariant(variant);

    std::stringstream attOut;
    for (const auto& att : gff.attributes_) {
        if (!attOut.str().empty()) attOut << ';';
        attOut << att.first << '=' << att.second;
    }

    std::stringstream line;
    line << gff.seqId_ << '\t' << gff.source_ << '\t' << gff.type_ << '\t' << gff.start_ << '\t'
         << gff.end_ << '\t' << gff.score_ << '\t' << gff.strand_ << '\t' << gff.phase_ << '\t'
         << attOut.str();
    WriteLine(line.str());
}

void GffWriter::WriteVariants(const std::vector<Variant>& variants)
{
    for (const auto& v : variants)
        WriteVariant(v);
}

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
