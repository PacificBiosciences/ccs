// Author: David Alexander

#include "BasicDenseMatrix.h"

#include <stdexcept>

namespace PacBio {
namespace Consensus {

BasicDenseMatrix::BasicDenseMatrix(size_t rows, size_t cols)
    : nCols_(cols), nRows_(rows), entries_(new double[nRows_ * nCols_])
{
}

BasicDenseMatrix::~BasicDenseMatrix() { delete[] entries_; }

void BasicDenseMatrix::ToHostMatrix(double **mat, int *rows, int *cols) const
{
    *mat = new double[Rows() * Columns()];
    *rows = Rows();
    *cols = Columns();
    for (size_t i = 0; i < Rows(); i++) {
        for (size_t j = 0; j < Columns(); j++) {
            (*mat)[i * Columns() + j] = (*this)(i, j);
        }
    }
}

size_t BasicDenseMatrix::UsedEntries() const { throw std::runtime_error("Unimplemented!"); }

float BasicDenseMatrix::UsedEntriesRatio() const { throw std::runtime_error("Unimplemented!"); }

size_t BasicDenseMatrix::AllocatedEntries() const { throw std::runtime_error("Unimplemented!"); }
}
}
