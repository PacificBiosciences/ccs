
#include "ScaledMatrix.h"

namespace PacBio {
namespace Consensus {

ScaledMatrix::ScaledMatrix(int rows, int cols) : SparseMatrix(rows, cols), logScalars_(cols, 0.0) {}
ScaledMatrix::ScaledMatrix(const ScaledMatrix& other)
    : SparseMatrix(other), logScalars_(other.logScalars_)
{
}

void ScaledMatrix::Reset(size_t rows, size_t cols)
{
    std::vector<double>(cols, 0.0).swap(logScalars_);
    SparseMatrix::Reset(rows, cols);
}

}  // namespace Consensus
}  // namespace PacBio
