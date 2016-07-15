
%{
#include <pacbio/consensus/Coverage.h>
%}

#ifdef SWIGPYTHON
    %apply (int DIM1, int* IN_ARRAY1)
         { (int tStartDim, int* tStart),
           (int tEndDim,   int* tEnd)  };
    %apply (int DIM1, int* ARGOUT_ARRAY1)
         { (int winLen, int* coverage) };
#endif // SWIGPYTHON

%include <pacbio/consensus/Coverage.h>
