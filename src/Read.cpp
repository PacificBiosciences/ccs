
#include <cassert>

#include <pacbio/consensus/ModelConfig.h>
#include <pacbio/consensus/Read.h>

namespace PacBio {
namespace Consensus {

SNR::SNR(const double a, const double c, const double g, const double t) : A(a), C(c), G(g), T(t) {}
SNR::SNR(const std::vector<double>& snrs) : A(snrs[0]), C(snrs[1]), G(snrs[2]), T(snrs[3])
{
    assert(snrs.size() == 4);
}

Read::Read(const std::string& name, const std::string& seq, const std::vector<uint8_t>& ipd,
           const std::vector<uint8_t>& pw, const SNR& snr, const std::string& model)
    : Name{name}, Seq{seq}, IPD{ipd}, PulseWidth{pw}, SignalToNoise{snr}, Model{model}
{
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
