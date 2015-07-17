%module ConsensusCore

%{
#define SWIG_FILE_WITH_INIT
%}

//
// SWIG doesn't support this even in languages with
// operator overloading.  Our vector-like objects
// need to expose an "ElementAt" accessor.
//
%ignore *::operator[];

%ignore boost::noncopyable;
namespace boost {
    class noncopyable {};
}

%include "Types.i"
#ifdef SWIGPYTHON
%include "numpy.i"
%numpy_typemaps(float, NPY_FLOAT, int)
#endif // SWIGPYTHON

%include "Evaluator.i"
%include "Exceptions.i"
%include "Features.i"
%include "Matrix.i"
%include "Mutation.i"
%include "QuiverConsensus.i"
%include "PairwiseAlignment.i"
%include "PoaConsensus.i"
%include "Utils.i"
%include "Statistics.i"
%include "Version.i"

%include "Edna.i"

%init %{
#ifdef SWIGPYTHON
  import_array();
#endif // SWIGPYTHON
%}
