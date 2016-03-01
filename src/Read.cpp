
#include <pacbio/consensus/Read.h>
#include <pacbio/consensus/ModelConfig.h>

namespace PacBio {
namespace Consensus {

Read::Read(const std::string& name, const std::string& seq, const std::vector<uint8_t>& cov,
           const std::string& model)
    : Name{name}, Seq{seq}, Cov{cov}, Model{model}
{
}

Read::Read(const std::string& name, const std::string& seq, const std::string& model)
    : Name{name}, Seq{seq}, Cov(seq.length(), 0), Model{model}
{
    for (size_t i = 0; i < seq.length(); ++i)
        Cov[i] = detail::TranslationTable[static_cast<unsigned char>(Seq[i])];
}

MappedRead::MappedRead(const Read& read, StrandEnum strand, size_t templateStart,
                       size_t templateEnd, bool pinStart, bool pinEnd)
    : Read(read)
    , Strand{strand}
    , TemplateStart{templateStart}
    , TemplateEnd{templateEnd}
    , PinStart{pinStart}
    , PinEnd{pinEnd}
{
}

std::ostream& operator<<(std::ostream& os, const MappedRead& mr)
{
    os << "MappedRead(Read(\"" << mr.Name << "\", \"" << mr.Seq << "\", \"" << mr.Model << "\"), ";
    if (mr.Strand == StrandEnum::FORWARD)
        os << "StrandEnum_FORWARD, ";
    else if (mr.Strand == StrandEnum::REVERSE)
        os << "StrandEnum_REVERSE, ";
    else if (mr.Strand == StrandEnum::UNMAPPED)
        os << "StrandEnum_UNMAPPED, ";
    os << mr.TemplateStart << ", " << mr.TemplateEnd << ", ";
    os << mr.PinStart << ", " << mr.PinEnd << ")";
    return os;
}

}  // namespace Consensus
}  // namespace PacBio
