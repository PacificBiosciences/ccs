// Author: Lance Hepler

#pragma once

#include <tuple>
#include <vector>

// Initialize data structures, do NOT remove
#include <pacbio/consensus/internal/ModelInternalInitializer.h>

#include <pacbio/consensus/Mutation.h>
#include <pacbio/consensus/PolishResult.h>

namespace PacBio {
namespace Consensus {

// forward declaration
class Integrator;

struct PolishConfig
{
    size_t MaximumIterations;
    size_t MutationSeparation;
    size_t MutationNeighborhood;

    bool Diploid;

    PolishConfig(size_t iterations = 40, size_t separation = 10, size_t neighborhood = 20,
                 bool diploid = false);
};

struct RepeatConfig
{
    size_t MaximumRepeatSize;
    size_t MinimumElementCount;
    size_t MaximumIterations;

    RepeatConfig(size_t repeatSize = 3, size_t elementCount = 3, size_t iterations = 40);
};

/// Given an Integrator and a PolishConfig,
/// iteratively polish the template,
/// and return meta information about the procedure.
///
/// The template will be polished within the Integrator.
PolishResult Polish(Integrator* ai, const PolishConfig& cfg);

PolishResult PolishRepeats(Integrator* ai, const RepeatConfig& cfg);

/// Struct that contains vectors for the base-wise individual and compound QVs.
struct QualityValues
{
    std::vector<int> Qualities;
    std::vector<int> DeletionQVs;
    std::vector<int> InsertionQVs;
    std::vector<int> SubstitutionQVs;
};

/// Generates phred qualities of the current template.
std::vector<int> ConsensusQualities(Integrator& ai);

/// Generates individual and compound phred qualities of the current template.
QualityValues ConsensusQVs(Integrator& ai);

/// Returns a list of all possible mutations that can be applied to the template
/// of the provided integrator.
std::vector<Mutation> Mutations(const Integrator& ai, bool diploid = false);

/// Returns a list of all possible repeat mutations of the template
/// of the provided integrator
std::vector<Mutation> RepeatMutations(const Integrator& ai, const RepeatConfig& cfg);

}  // namespace Consensus
}  // namespace PacBio
