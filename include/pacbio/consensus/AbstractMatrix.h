// Author: David Alexander

#pragma once

#include <cstddef>

// Initialize data structures, do NOT remove
#include <pacbio/consensus/internal/ModelInternalInitializer.h>

namespace PacBio {
namespace Consensus {

/// AbstractMatrix is a superclass of the matrix types used in the arrow
/// banded dynamic programming.  It exposes a minimal interface only
/// intended for diagnostic purposes (looking at a matrix from Python,
/// seeing how well the banding is working, ...).  No matrix implementation
/// details are exposed---one can think of this as effectively an opaque
/// data type.
class AbstractMatrix
{
public:
    /// Method SWIG clients can use to get a native matrix (e.g. Numpy)
    /// mat must be filled as a ROW major matrix.
    /// N.B.: Needs int, not size_t dimensions, for SWIG/numpy
    virtual void ToHostMatrix(double** mat, int* rows, int* cols) const = 0;

public:
    // Methods for inquiring about matrix occupancy.
    virtual size_t UsedEntries() const = 0;
    virtual float UsedEntriesRatio() const = 0;
    virtual size_t AllocatedEntries() const = 0;  // an entry may be allocated but not used

public:
    virtual ~AbstractMatrix() {}
};

}  // namespace Consensus
}  // namespace PacBio
