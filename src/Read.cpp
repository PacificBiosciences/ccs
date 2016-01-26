
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

}  // namespace Consensus
}  // namespace PacBio
