// Author: Lance Hepler

#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

#include "SparseMatrix.h"

namespace PacBio {
namespace Consensus {

/// This class inherits from SparseMatrix and extends it by having a
/// column-wise scaling factor.
class ScaledMatrix : public SparseMatrix
{
public:
    enum Direction
    {
        FORWARD,
        REVERSE
    };

public:  // constructor/destructor
    /// Constructor with explicit dimensions.
    ScaledMatrix(size_t rows, size_t cols, Direction dir);
    /// Copy constructor.
    ScaledMatrix(const ScaledMatrix& other);
    /// Destructor.
    ~ScaledMatrix(void) override = default;

public:
    /// Clears and resizes the internal data structures.
    void Reset(size_t rows, size_t cols) override;
    /// Set direction and reset column-wise log scalars.
    Direction SetDirection(Direction dir);

public:  // nullability
    /// Returns a ScaledMatrix representing null.
    static const ScaledMatrix& Null();

public:  // information about entries filled by column
    /// Rescale column j by max_val.
    /// If maxProvided is false, determine max_val.
    template <bool maxProvided>
    void FinishEditingColumn(size_t j, size_t usedBegin, size_t usedEnd, double max_val = 0.0);

public:  // Scaling and normalization
    /// Get the log scale for column j.
    double GetLogScale(size_t j) const;
    /// Get the log scales from colum s to e
    double GetLogProdScales(size_t s, size_t e) const;
    /// Get the log scale
    double GetLogProdScales() const;

public:  // Convenient matrix access for SWIG
    /// Convert sparse to full matrix.
    void ToHostMatrix(double** mat, int* rows, int* cols) const override;

private:
    std::vector<double> logScalars_;
    Direction dir_;
};

std::ostream& operator<<(std::ostream&, const ScaledMatrix&);

inline const ScaledMatrix& ScaledMatrix::Null()
{
    static auto nullObj = std::make_unique<ScaledMatrix>(0, 0, FORWARD);
    return *nullObj;
}

template <bool maxProvided>
inline void ScaledMatrix::FinishEditingColumn(const size_t j, const size_t usedBegin,
                                              const size_t usedEnd, double max_val)
{
    // get the constant to scale by
    if (!maxProvided) {
        max_val = 0.0;
        for (size_t i = usedBegin; i < usedEnd; ++i) {
            max_val = std::max(max_val, SparseMatrix::Get(i, j));
        }
    }

    // cumsum stuff
    double last = 0.0;
    if (dir_ == FORWARD && j > 0)
        last = logScalars_[j - 1];
    else if (dir_ == REVERSE && j + 1 < Columns())
        last = logScalars_[j + 1];

    // set it
    if (max_val != 0.0 && max_val != 1.0) {
        for (size_t i = usedBegin; i < usedEnd; ++i) {
            SparseMatrix::Set(i, j, SparseMatrix::Get(i, j) / max_val);
        }
        logScalars_[j] = last + std::log(max_val);
    } else {
        logScalars_[j] = last;
    }

    SparseMatrix::FinishEditingColumn(j, usedBegin, usedEnd);
}

inline double ScaledMatrix::GetLogScale(size_t j) const { return logScalars_[j]; }

inline double ScaledMatrix::GetLogProdScales(size_t beginColumn, size_t endColumn) const
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
