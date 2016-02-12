
#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

#include "SparseMatrix.h"

namespace PacBio {
namespace Consensus {

class ScaledMatrix : public SparseMatrix
{
public:
    enum Direction
    {
        FORWARD,
        REVERSE
    };

public:  // constructor/destructor
    ScaledMatrix(int rows, int cols, Direction dir);
    ScaledMatrix(const ScaledMatrix& other);
    ~ScaledMatrix(void) = default;

public:
    void Reset(size_t rows, size_t cols);
    Direction SetDirection(Direction dir);

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
    Direction dir_;
};

std::ostream& operator<<(std::ostream&, const ScaledMatrix&);

inline const ScaledMatrix& ScaledMatrix::Null()
{
    static ScaledMatrix* nullObj = new ScaledMatrix(0, 0, FORWARD);
    return *nullObj;
}

inline void ScaledMatrix::FinishEditingColumn(const int j, const int usedBegin, const int usedEnd)
{
    // get the constant to scale by
    double c = 0.0;
    for (int i = usedBegin; i < usedEnd; ++i) {
        c = std::max(c, SparseMatrix::Get(i, j));
    }

    // cumsum stuff
    double last = 0.0;
    if (dir_ == FORWARD && j > 0)
        last = logScalars_[j - 1];
    else if (dir_ == REVERSE && j + 1 < Columns())
        last = logScalars_[j + 1];

    // set it
    if (c != 0.0 && c != 1.0) {
        for (int i = usedBegin; i < usedEnd; ++i) {
            SparseMatrix::Set(i, j, SparseMatrix::Get(i, j) / c);
        }
        logScalars_[j] = last + std::log(c);
    } else {
        logScalars_[j] = last;
    }

    SparseMatrix::FinishEditingColumn(j, usedBegin, usedEnd);
}

inline double ScaledMatrix::GetLogScale(int j) const { return logScalars_[j]; }
inline double ScaledMatrix::GetLogProdScales(int beginColumn, int endColumn) const
{
    double f, l;
    if (dir_ == FORWARD) {
        f = (beginColumn > 0) ? logScalars_[beginColumn - 1] : 0.0;
        l = (endColumn > 0) ? logScalars_[endColumn - 1] : 0.0;
        return l - f;
    }  // else dir_ == REVERSE
    f = logScalars_[beginColumn];
    l = (endColumn < Columns()) ? logScalars_[endColumn] : 0.0;
    return f - l;
}

inline double ScaledMatrix::GetLogProdScales() const
{
    if (dir_ == FORWARD) return logScalars_.back();
    return logScalars_.front();
}

}  // namespace Consensus
}  // namespace PacBio
