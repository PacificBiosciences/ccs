// Author: David Alexander

#pragma once

// Initialize data structures, do NOT remove
#include <pacbio/consensus/internal/ModelInternalInitializer.h>

namespace PacBio {

namespace Data {
struct MappedRead;
}  // PacBio::Data

namespace Consensus {

class Evaluator;

//
// EasyReadScorer is a simplified means to get access to the evaluator
// scoring machinery.  This is suitable for experimental purposes,
// not intended for use in production code.
//
class EasyReadScorer
{

public:
    //
    // This is a convenient means to build an evaluator---it will automatically
    // look up the model, populate the template object, etc.
    //
    static Evaluator* MakeEvaluator(const std::string& tplString,
                                    const PacBio::Data::MappedRead& mappedRead, double minZScore,
                                    double scoreDiff);
};

}  // PacBio::Consensus
}  // PacBio
