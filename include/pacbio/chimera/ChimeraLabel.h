// Author: Armin Töpfer

#pragma once

#include <stdbool.h>
#include <string>

namespace PacBio {
namespace Chimera {

///
/// \brief Label that annotates a read for a single chimeric breakpoint
///
struct ChimeraLabel
{
    // Instance variables
    std::string sequenceId;
    bool chimeraFlag;
    std::string leftParentId;
    std::string rightParentId;
    int32_t crossover;
    double score;

    // Default Constructor
    ChimeraLabel(std::string sequenceIdArg, std::string leftParentArg, std::string rightParentArg,
                 int32_t crossoverArg, double scoreArg)
        : sequenceId(sequenceIdArg)
        , chimeraFlag(false)
        , leftParentId(leftParentArg)
        , rightParentId(rightParentArg)
        , crossover(crossoverArg)
        , score(scoreArg){};

    // Name-Only or Place-Holder Constructor
    explicit ChimeraLabel(std::string sequenceIdArg)
        : sequenceId(sequenceIdArg)
        , chimeraFlag(false)
        , leftParentId("N/A")
        , rightParentId("N/A")
        , crossover(-1)
        , score(-1.0){};

    // Empty or Dummy Constructor
    ChimeraLabel()
        : sequenceId("Dummy")
        , chimeraFlag(false)
        , leftParentId("N/A")
        , rightParentId("N/A")
        , crossover(-1)
        , score(-1.0){};

    // Move constructor
    ChimeraLabel(ChimeraLabel &&src) = default;
    // Copy constructor is deleted!
    ChimeraLabel(const ChimeraLabel &src) = default;
    // Move assignment constructor
    ChimeraLabel &operator=(ChimeraLabel &&rhs) = default;
    // Copy assignment constructor is deleted!
    ChimeraLabel &operator=(const ChimeraLabel &rhs) = default;
    // Destructor
    ~ChimeraLabel() = default;

    friend std::ostream &operator<<(std::ostream &o, const ChimeraLabel &label)
    {
        // Stream the Sequence Id first
        o << label.sequenceId << ",";

        // Then a human-readable representation of the flag
        if (label.chimeraFlag)
            o << "True"
              << ",";
        else
            o << "False"
              << ",";

        // The score is only meaningfully defined > 0
        if (label.score >= 0.0)
            o << label.score << ",";
        else
            o << "NaN"
              << ",";

        // Finally the parents and the putative crossover
        o << label.leftParentId << "," << label.rightParentId << "," << label.crossover;

        // Return the stream reference
        return o;
    }
};

}  // namespace Chimera
}  // namespace PacBio
