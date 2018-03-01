// Author: Lance Hepler

#pragma once

#include <array>
#include <cstdint>
#include <map>
#include <memory>
#include <ostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include <pacbio/UnanimityConfig.h>

// Initialize data structures, do NOT remove
#include <pacbio/consensus/internal/ModelInternalInitializer.h>

#include <pacbio/data/Read.h>
#include <pacbio/data/internal/BaseEncoding.h>

namespace PacBio {

namespace Data {
struct MappedRead;
struct SNR;
}

namespace Consensus {

// fwd decl
class AbstractRecursor;
class AbstractTemplate;

using Data::detail::NCBI2na;
using Data::detail::NCBI4na;

// The allele representation in Unanimity
// currently employs the NCBI4na model, in
// order to account for diploid sites.
typedef NCBI4na AlleleRep;

struct TemplatePosition
{
    char Base;
    AlleleRep Idx;
    double Match;
    double Branch;
    double Stick;
    double Deletion;

public:
    // Constructor for backwards-compatibility
    UNANIMITY_CONSTEXPR TemplatePosition(const char base, const double match, const double branch,
                                         const double stick, const double deletion)
        : Base{base}
        , Idx{AlleleRep::FromASCII(Base)}
        , Match{match}
        , Branch{branch}
        , Stick{stick}
        , Deletion{deletion}
    {
    }
};

std::ostream& operator<<(std::ostream&, const TemplatePosition&);

enum struct MoveType : uint8_t
{
    MATCH = 0,
    BRANCH = 1,
    STICK = 2,
    DELETION = 3  // never used for covariate
};

enum struct MomentType : uint8_t
{
    FIRST = 0,
    SECOND = 1
};

class ModelConfig
{
public:
    virtual ~ModelConfig() {}
    virtual std::unique_ptr<AbstractRecursor> CreateRecursor(const PacBio::Data::MappedRead& mr,
                                                             double scoreDiff) const = 0;
    virtual std::vector<TemplatePosition> Populate(const std::string& tpl) const = 0;
    virtual std::pair<Data::Read, std::vector<MoveType>> SimulateRead(
        std::default_random_engine* const rng, const std::string& tpl,
        const std::string& readname) const = 0;
    virtual double ExpectedLLForEmission(MoveType move, const AlleleRep& prev,
                                         const AlleleRep& curr, MomentType moment) const = 0;
};

}  // namespace Consensus
}  // namespace PacBio
