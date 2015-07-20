%{
/* Includes the header in the wrapper code */
#include <ConsensusCore/Feature.hpp>
#include <ConsensusCore/Features.hpp>
using namespace ConsensusCore;
%}

#if SWIGPYTHON

%include "numpy.i"
%numpy_typemaps(float, NPY_FLOAT, int)

%apply (float* IN_ARRAY1, int DIM1)
     { (const float* inPtr, int length) };

#endif // SWIGPYTHON

#if SWIGCSHARP
%include "arrays_csharp.i"
%apply unsigned char INPUT[] { unsigned char* }
%csmethodmodifiers ConsensusCore::Feature::ToString() const "public override"
#endif // SWIGCSHARP


%include <ConsensusCore/Feature.hpp>
%include <ConsensusCore/Features.hpp>

namespace ConsensusCore {
    %template(CharFeature) Feature<char>;
    %template(FloatFeature) Feature<float>;
    %template(IntFeature) Feature<int>;
}

%include "carrays.i"
%array_class(float, FloatArray);
%array_class(int, IntArray);
