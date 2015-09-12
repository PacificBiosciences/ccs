
#pragma once

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

#include "SparseMatrix.h"

namespace PacBio {
namespace Consensus {

class ScaledMatrix : public SparseMatrix
{
public:  // constructor/destructor
    ScaledMatrix(int rows, int cols);
    ScaledMatrix(const ScaledMatrix& other);
    ~ScaledMatrix(void) = default;

public:
    void Reset(size_t rows, size_t cols);

public:  // nullability
    static const ScaledMatrix& Null();

public:  // information about entries filled by column
    void FinishEditingColumn(int j, int usedBegin, int usedEnd);

public:  // Scaling and normalization
    double GetLogScale(int j) const;
    double GetLogProdScales(int s, int e) const;
    double GetLogProdScales() const;

private:
    std::vector<double> logScalars_;
};

inline
const ScaledMatrix& ScaledMatrix::Null()
{
    static ScaledMatrix* nullObj = new ScaledMatrix(0, 0);
    return *nullObj;
}

inline
void ScaledMatrix::FinishEditingColumn(int j, int usedBegin, int usedEnd)
{
    // get the constant to scale by
    double c = 0.0;
    for (int i = usedBegin; i < usedEnd; ++i)
    {
        c = std::max(c, SparseMatrix::Get(i, j));
    }

    // set it
    if (c != 0.0 && c != 1.0)
    {
        for (int i = usedBegin; i < usedEnd; ++i)
        {
            SparseMatrix::Set(i, j, SparseMatrix::Get(i, j) / c);
        }
        logScalars_[j] = std::log(c);
    }
    else
    {
        logScalars_[j] = 0.0;
    }

    SparseMatrix::FinishEditingColumn(j, usedBegin, usedEnd);
}

inline
double ScaledMatrix::GetLogScale(int j) const
{
    return logScalars_[j];
}

inline
double ScaledMatrix::GetLogProdScales(int beginColumn, int endColumn) const
{
    return std::accumulate(logScalars_.begin() + beginColumn,
                           logScalars_.begin() + endColumn,
                           0.0);
}

inline
double ScaledMatrix::GetLogProdScales() const
{
    return std::accumulate(logScalars_.begin(),
                           logScalars_.end(),
                           0.0);
}

} // namespace Consensus
} // namespace PacBio
