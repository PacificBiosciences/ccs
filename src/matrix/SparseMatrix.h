// Author: David Alexander

#pragma once

#include <algorithm>
#include <cassert>
#include <memory>
#include <utility>
#include <vector>

#include <pacbio/consensus/AbstractMatrix.h>
#include "SparseVector.h"

namespace PacBio {
namespace Consensus {

/// The SparseMatrix is based on a vector of SparseVectors.
class SparseMatrix : public AbstractMatrix
{
public:  // Constructor, destructor
    /// Constructor with explicit dimensions.
    SparseMatrix(size_t rows, size_t cols);
    /// Copy constructor.
    SparseMatrix(const SparseMatrix& other);
    /// Destructor.
    virtual ~SparseMatrix();

public:
    /// Clears and resizes the internal data structures.
    virtual void Reset(size_t rows, size_t cols);

public:  // Nullability
    /// Returns a SparseMatrix representing null.
    static const SparseMatrix& Null();
    /// Returns if the both dimensions are zero.
    bool IsNull() const;

public:  // Size information
    size_t Rows() const;
    size_t Columns() const;

public:  // Information about entries filled by column
    /// Prepared the underlying SparseVector of column j
    /// from rows hintBegin to hintEnd.
    void StartEditingColumn(size_t j, size_t hintBegin, size_t hintEnd);
    /// Finish editing column j and store the used rows in usedRanges_.
    void FinishEditingColumn(size_t j, size_t usedBegin, size_t usedEnd);
    /// Retreive the row range for column j.
    std::pair<size_t, size_t> UsedRowRange(size_t j) const;
    // Checks if no rows are set for column j.
    bool IsColumnEmpty(size_t j) const;
    /// Computes the number of filled cells.
    size_t UsedEntries() const override;
    /// Computes the ratio of filled cells.
    float UsedEntriesRatio() const override;
    /// Computes the number of allocated cells.
    /// An entry may be allocated but not used.
    size_t AllocatedEntries() const override;

public:  // Accessors
    /// Access cell at row i and column j.
    /// If not allocated, return 0.
    const double& operator()(size_t i, size_t j) const;
    /// Checks if cell is allocated.
    bool IsAllocated(size_t i, size_t j) const;
    double Get(size_t i, size_t j) const;
    void Set(size_t i, size_t j, double v);
    /// Clear content of column j and reset respective row range.
    void ClearColumn(size_t j);

public:
    /// Convert sparse to full matrix.
    void ToHostMatrix(double** mat, int* rows, int* cols) const override;

private:
    static double* EmptyCell();

private:
    void CheckInvariants(size_t column) const;

private:
    std::vector<std::unique_ptr<SparseVector>> columns_;
    size_t nCols_;
    size_t nRows_;
    size_t columnBeingEdited_;
    std::vector<std::pair<size_t, size_t>> usedRanges_;
};

//
// Nullability
//
inline const SparseMatrix& SparseMatrix::Null()
{
    static auto nullObj = std::make_unique<SparseMatrix>(0, 0);
    return *nullObj;
}

inline bool SparseMatrix::IsNull() const { return (Rows() == 0 && Columns() == 0); }

//
// Size information
//
inline size_t SparseMatrix::Rows() const { return nRows_; }
inline size_t SparseMatrix::Columns() const { return nCols_; }

//
// Entry range queries per column
//
inline void SparseMatrix::StartEditingColumn(size_t j, size_t hintBegin, size_t hintEnd)
{
    assert(columnBeingEdited_ == std::numeric_limits<size_t>::max());
    columnBeingEdited_ = j;
    if (columns_[j] != NULL) {
        columns_[j]->ResetForRange(hintBegin, hintEnd);
    } else {
        columns_[j] = std::make_unique<SparseVector>(Rows(), hintBegin, hintEnd);
    }
}

inline void SparseMatrix::FinishEditingColumn(size_t j, size_t usedRowsBegin, size_t usedRowsEnd)
{
    assert(columnBeingEdited_ == j);
    usedRanges_[j] = std::make_pair(usedRowsBegin, usedRowsEnd);
    CheckInvariants(columnBeingEdited_);
    columnBeingEdited_ = std::numeric_limits<size_t>::max();
}

inline std::pair<size_t, size_t> SparseMatrix::UsedRowRange(size_t j) const
{
    assert(j < usedRanges_.size());
    return usedRanges_[j];
}

inline bool SparseMatrix::IsColumnEmpty(size_t j) const
{
    assert(j < usedRanges_.size());
    size_t begin, end;
    std::tie(begin, end) = usedRanges_[j];
    return begin >= end;
}

//
// Accessors
//
inline const double& SparseMatrix::operator()(size_t i, size_t j) const
{
    if (columns_[j] == NULL) {
        static const double emptyCell = 0.0;
        return emptyCell;
    } else {
        return (*columns_[j])(i);
    }
}

inline bool SparseMatrix::IsAllocated(size_t i, size_t j) const
{
    return columns_[j] != NULL && columns_[j]->IsAllocated(i);
}

inline double SparseMatrix::Get(size_t i, size_t j) const { return (*this)(i, j); }

inline void SparseMatrix::Set(size_t i, size_t j, double v)
{
    assert(columnBeingEdited_ == j);
    columns_[j]->Set(i, v);
}

inline void SparseMatrix::ClearColumn(size_t j)
{
    usedRanges_[j] = std::make_pair(0, 0);
    columns_[j]->Clear();
    CheckInvariants(j);
}

}  // namespace Consensus
}  // namespace PacBio
