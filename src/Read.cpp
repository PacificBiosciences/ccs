
#include <pacbio/consensus/Read.h>

namespace PacBio {
namespace Consensus {

Read::Read(const std::string& name,
           const std::string& seq,
#ifndef NO_COVARIATES
           const std::vector<uint8_t>& cov,
#endif // NO_COVARIATES
           const std::string& model)
    : name_{name}
    , seq_{seq}
#ifndef NO_COVARIATES
    , cov_{cov}
#endif // NO_COVARIATES
    , mdl_{model}
{ }

const std::string& Read::Name() const
{
    return name_;
}

const std::string& Read::Seq() const
{
    return seq_;
}

#ifndef NO_COVARIATES
const std::vector<uint8_t>& Read::Cov() const
{
    return cov_;
}
#endif // NO_COVARIATES

const std::string& Read::Model() const
{
    return mdl_;
}

MappedRead::MappedRead(const Read& read,
                       StrandEnum strand,
                       size_t templateStart,
                       size_t templateEnd,
                       bool pinStart,
                       bool pinEnd)
    : Read(read)
    , strand_{strand}
    , tplStart_{templateStart}
    , tplEnd_{templateEnd}
    , pinStart_{pinStart}
    , pinEnd_{pinEnd}
{ }

StrandEnum MappedRead::Strand() const
{
    return strand_;
}

size_t MappedRead::TemplateStart() const
{
    return tplStart_;
}

size_t MappedRead::TemplateEnd() const
{
    return tplEnd_;
}

bool MappedRead::PinStart() const
{
    return pinStart_;
}

bool MappedRead::PinEnd() const
{
    return pinEnd_;
}

} // namespace Consensus
} // namespace PacBio
