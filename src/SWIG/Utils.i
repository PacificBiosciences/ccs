%{
/* Includes the header in the wrapper code */
#include <ConsensusCore/Utils.hpp>
#include <ConsensusCore/Coverage.hpp>
#include <ConsensusCore/Logging.hpp>
using namespace ConsensusCore;
%}

%include <ConsensusCore/Types.hpp>

#ifdef SWIGPYTHON
    %apply (int DIM1, int* IN_ARRAY1)
         { (int tStartDim, int* tStart),
           (int tEndDim,   int* tEnd)  };
    %apply (int DIM1, int* ARGOUT_ARRAY1)
         { (int winLen, int* coverage) };
#endif // SWIGPYTHON

%include <ConsensusCore/Utils.hpp>
%include <ConsensusCore/Coverage.hpp>
%include <ConsensusCore/Logging.hpp>
