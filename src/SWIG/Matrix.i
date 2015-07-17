%{
/* Includes the header in the wrapper code */
#include <ConsensusCore/Types.hpp>
#include <ConsensusCore/Matrix/AbstractMatrix.hpp>
#include <ConsensusCore/Matrix/DenseMatrix.hpp>
#include <ConsensusCore/Matrix/SparseMatrix.hpp>
using namespace ConsensusCore;
%}

%include <ConsensusCore/Types.hpp>

#ifdef SWIGPYTHON
        // apply this typemap to ToHostMatrix
        %apply (float** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2)
             { (float** mat, int* rows, int* cols) };
#endif // SWIGPYTHON

%newobject *::UsedRowRange;


%include <ConsensusCore/Matrix/AbstractMatrix.hpp>
%include <ConsensusCore/Matrix/DenseMatrix.hpp>
%include <ConsensusCore/Matrix/SparseMatrix.hpp>
