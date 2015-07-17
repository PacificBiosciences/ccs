
#pragma once

#include <vector>

#include <ConsensusCore/Matrix/SparseMatrix.hpp>

namespace ConsensusCore {

    template<typename M>
    class ScaledMatrix : public M
    {
    private:
        typedef typename M::FloatType F;

    public:  // constructor/destructor
        ScaledMatrix(int rows, int cols);
        ScaledMatrix(const ScaledMatrix& other);
        ~ScaledMatrix(void) = default;

    public:  // nullability
        static const ScaledMatrix<M>& Null();

    public:  // information about entries filled by column
        void FinishEditingColumn(int j, int usedBegin, int usedEnd);

    public:  // Scaling and normalization
        F GetLogScale(int j) const;
        F GetLogProdScales(int s, int e) const;
        F GetLogProdScales() const;

    private:
        std::vector<F> logScalars_;
    };

    typedef ScaledMatrix<SparseMatrix<double, double>> ScaledSparseMatrixD;
}

#include "ScaledMatrix-inl.hpp"
