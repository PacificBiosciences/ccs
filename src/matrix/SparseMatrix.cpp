// Author: David Alexander

#include "SparseMatrix.h"

#include <memory>

namespace PacBio {
namespace Consensus {

SparseMatrix::SparseMatrix(const size_t rows, const size_t cols)
    : columns_(cols)
    , nCols_(cols)
    , nRows_(rows)
    , columnBeingEdited_(std::numeric_limits<size_t>::max())
    , usedRanges_(cols, std::make_pair(0, 0))
{
}

SparseMatrix::SparseMatrix(const SparseMatrix& other)
    : columns_(other.nCols_)
    , nCols_(other.nCols_)
    , nRows_(other.nRows_)
    , columnBeingEdited_(other.columnBeingEdited_)
    , usedRanges_(other.usedRanges_)
{
    for (size_t j = 0; j < nCols_; j++)
        if (other.columns_[j]) columns_[j] = std::make_unique<SparseVector>(*other.columns_[j]);
}

SparseMatrix::~SparseMatrix() = default;

void SparseMatrix::Reset(const size_t rows, const size_t cols)
{
    std::vector<std::unique_ptr<SparseVector>>(cols).swap(columns_);
    nCols_ = cols;
    nRows_ = rows;
    std::vector<std::pair<size_t, size_t>>(cols, std::make_pair(0, 0)).swap(usedRanges_);
    columnBeingEdited_ = std::numeric_limits<size_t>::max();
}

size_t SparseMatrix::UsedEntries() const
{
    // use column ranges
    size_t filledEntries = 0;
    for (size_t col = 0; col < Columns(); ++col) {
        size_t start, end;
        std::tie(start, end) = UsedRowRange(col);
        filledEntries += (end - start);
    }
    return filledEntries;
}

float SparseMatrix::UsedEntriesRatio() const
{
    const auto filled = static_cast<float>(UsedEntries());
    const auto size = static_cast<float>(Rows() * Columns());
    return filled / size;
}

size_t SparseMatrix::AllocatedEntries() const
{
    size_t sum = 0;
    for (size_t j = 0; j < nCols_; j++) {
        sum += (columns_[j] != nullptr ? columns_[j]->AllocatedEntries() : 0);
    }
    return sum;
}

void SparseMatrix::ToHostMatrix(double** mat, int* rows, int* cols) const
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    *mat = new double[Rows() * Columns()];
    *rows = static_cast<int>(Rows());
    *cols = static_cast<int>(Columns());
    for (size_t i = 0; i < Rows(); i++) {
        for (size_t j = 0; j < Columns(); j++) {
            (*mat)[i * Columns() + j] = IsAllocated(i, j) ? Get(i, j) : nan;
        }
    }
}

void SparseMatrix::CheckInvariants(size_t column) const
{
#ifndef NDEBUG
    for (size_t j = 0; j < nCols_; j++) {
        if (columns_[j] != nullptr) columns_[j]->CheckInvariants();
    }
#endif
}

}  // namespace Consensus
}  // namespace PacBio
