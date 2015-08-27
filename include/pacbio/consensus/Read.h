
#pragma once

namespace PacBio {
namespace Consensus {

class Read
{
public:
    Read(const std::string& name,
         const std::string& seq,
#ifndef NO_COVARIATES
         const std::vector<uint8_t>& cov,
#endif // NO_COVARIATES
         const std::string& model);

    // TODO(lhepler): inlines?
    const std::string& Name() const;
    const std::string& Seq() const;
#ifndef NO_COVARIATES
    const std::vector<uint8_t>& Cov() const;
#endif // NO_COVARIATES
    const std::string& Model() const;

private:
    std::string name_;
    std::string seq_;
#ifndef NO_COVARIATES
    std::vector<uint8_t> cov_;
#endif // NO_COVARIATES
    std::string mdl_;
};

enum StrandEnum
{
    FORWARD,
    REVERSE
};

class MappedRead : public Read
{
public:
    MappedRead(const Read& read,
               StrandEnum strand,
               size_t templateStart,
               size_t templateEnd,
               bool pinStart = false,
               bool pinEnd = false);

    // TODO(lhepler): inlines?
    StrandEnum Strand() const;
    size_t TemplateStart() const;
    size_t TemplateEnd() const;
    bool PinStart() const;
    bool PinEnd() const;

private:
    StrandEnum strand_;
    size_t tplStart_;
    size_t tplEnd_;
    bool pinStart_;
    bool pinEnd_;
};

} // namespace Consensus
} // namespace PacBio
