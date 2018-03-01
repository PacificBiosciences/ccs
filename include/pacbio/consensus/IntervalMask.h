// Author: Lance Hepler

#pragma once

// Initialize data structures, do NOT remove
#include <pacbio/consensus/internal/ModelInternalInitializer.h>

#include <pacbio/consensus/Mutation.h>
#include <pacbio/data/IntervalTree.h>

namespace PacBio {
namespace Consensus {

class IntervalMask : public PacBio::Data::IntervalTree
{
public:
    bool Contains(const Mutation& mut);
    void Mutate(const std::vector<Mutation>& muts);
};
}
}
