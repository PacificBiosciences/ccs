// Author: Armin Töpfer

#pragma once

#include <vector>

#include <cstddef>

#include <boost/optional.hpp>

// Initialize data structures, do NOT remove
#include <pacbio/consensus/internal/ModelInternalInitializer.h>

#include <pacbio/consensus/Mutation.h>

namespace PacBio {
namespace Consensus {

struct DiploidSite
{
    MutationType mutType;
    std::vector<char> mutants;
    int64_t pos;
    boost::optional<double> pvalue;

    DiploidSite(const MutationType mutType_, const std::vector<char>& mutants_, const int64_t pos_,
                const boost::optional<double> pvalue_ = boost::none)
        : mutType{mutType_}, mutants{mutants_}, pos{pos_}, pvalue{pvalue_}
    {
    }

    bool operator==(const DiploidSite& rhs) const
    {
        return std::tie(mutType, mutants, pos, pvalue) ==
               std::tie(rhs.mutType, rhs.mutants, rhs.pos, rhs.pvalue);
    }
};

/// This struct contains the results of Integrator::Polish()
struct PolishResult
{
    // Did Polish() converge?
    bool hasConverged = false;
    // How many mutations have been tested?
    size_t mutationsTested = 0;
    // How many mutations have been actually applied?
    size_t mutationsApplied = 0;

    // For each iteration in Polish(), get the max of all Evaluators to
    // diagnose the worst performing one.
    //
    // Maximal ratio of populated alpha cells
    std::vector<float> maxAlphaPopulated;
    // Maximal ratio of populated beta cells
    std::vector<float> maxBetaPopulated;
    // Maximal number of flip flop events
    std::vector<int> maxNumFlipFlops;

    // Diploid results
    // The vector is sorted according to the standard
    // unanimity Mutation class criterion
    std::vector<DiploidSite> diploidSites;
};

PolishResult operator+(const PolishResult& lhs, const PolishResult& rhs);
}
}  // ::PacBio::Consensus