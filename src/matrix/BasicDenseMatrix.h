// Author: David Alexander

#pragma once

#include <cstddef>
#include "pacbio/consensus/AbstractMatrix.h"

namespace PacBio {
namespace Consensus {

//
// BasicDenseMatrix is a *basic* dense matrix, for use as an
// intermediate in matrix viewing operations (not in production code).
//
// It does not fully implement the interface that would be required to
// drop it in as a replacement for ScaledMatrix in the production code
// (i.e. for the recursor, etc.).  ConsensusCore did offer such a
// matrix, "DenseMatrix", and we could consider resurrecting such a
// class in the future.
//
class BasicDenseMatrix : public AbstractMatrix
{
private:
    size_t nCols_;
    size_t nRows_;
    double* entries_;

public:  // structors
    BasicDenseMatrix(size_t rows, size_t cols);
    BasicDenseMatrix(const BasicDenseMatrix& other) = delete;
    ~BasicDenseMatrix();

public:  // Size information
    size_t Rows() const;
    size_t Columns() const;

public:  // accessors
    double& operator()(size_t i, size_t j) const;

public:  // AbstractMatrix interface
    void ToHostMatrix(double** mat, int* rows, int* cols) const;
    size_t UsedEntries() const;
    float UsedEntriesRatio() const;
    size_t AllocatedEntries() const;
};

inline size_t BasicDenseMatrix::Rows() const { return nRows_; }

inline size_t BasicDenseMatrix::Columns() const { return nCols_; }

inline double& BasicDenseMatrix::operator()(size_t i, size_t j) const
{
    return entries_[i * Columns() + j];
}
}
}
