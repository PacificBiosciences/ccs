// Author: Lance Hepler

#include <cassert>
#include <utility>

#include <pacbio/consensus/ModelConfig.h>
#include <pacbio/data/Read.h>

namespace PacBio {
namespace Data {

SNR::SNR(const double a, const double c, const double g, const double t) : A(a), C(c), G(g), T(t) {}

SNR::SNR(const std::vector<float>& snrs) : A(snrs[0]), C(snrs[1]), G(snrs[2]), T(snrs[3])
{
    assert(snrs.size() == 4);
}

SNR::SNR(const std::vector<double>& snrs) : A(snrs[0]), C(snrs[1]), G(snrs[2]), T(snrs[3])
{
    assert(snrs.size() == 4);
}

namespace {

double clamp(double val, double lo, double hi) { return std::min(std::max(val, lo), hi); }

}  // namespace

SNR ClampSNR(const SNR& val, const SNR& lo, const SNR& hi)
{
    return SNR(clamp(val.A, lo.A, hi.A), clamp(val.C, lo.C, hi.C), clamp(val.G, lo.G, hi.G),
               clamp(val.T, lo.T, hi.T));
}

Read::Read(const std::string& name, const std::string& seq, const std::vector<uint8_t>& ipd,
           const std::vector<uint8_t>& pw, const SNR& snr, std::string model)
    : Name{name}, Seq{seq}, IPD{ipd}, PulseWidth{pw}, SignalToNoise{snr}, Model{std::move(model)}
{
    if (ipd.size() != seq.size() || pw.size() != seq.size()) {
        throw std::invalid_argument("Invalid Read (name=" + name +
                                    "): features IPD/PW/seq are of mismatched length");
    }
}

MappedRead::MappedRead(const Read& read, StrandType strand, size_t templateStart,
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
    switch (mr.Strand) {
        case StrandType::FORWARD:
            os << "StrandType_FORWARD, ";
            break;
        case StrandType::REVERSE:
            os << "StrandType_REVERSE, ";
            break;
        case StrandType::UNMAPPED:
            os << "StrandType_UNMAPPED, ";
            break;
        default:
            throw std::runtime_error("Unsupported Strand");
    }
    os << mr.TemplateStart << ", " << mr.TemplateEnd << ", ";
    os << mr.PinStart << ", " << mr.PinEnd << ")";
    return os;
}

}  // namespace Data
}  // namespace PacBio
