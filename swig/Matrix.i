%{
/* Includes the header in the wrapper code */
#include <pacbio/consensus/AbstractMatrix.h>
#include <pacbio/consensus/MatrixViewConvention.h>
%}

#ifdef SWIGPYTHON
        // apply this typemap to ToHostMatrix
        %apply (double** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2)
             { (double** mat, int* rows, int* cols) };
#endif // SWIGPYTHON

%include <pacbio/consensus/AbstractMatrix.h>
%include <pacbio/consensus/MatrixViewConvention.h>
