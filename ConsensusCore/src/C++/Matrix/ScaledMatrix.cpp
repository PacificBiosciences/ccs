
#include <ConsensusCore/Matrix/ScaledMatrix.hpp>
#include <ConsensusCore/Matrix/SparseMatrix.hpp>

namespace ConsensusCore {
    template class ScaledMatrix<SparseMatrix<double, double>>;
}
