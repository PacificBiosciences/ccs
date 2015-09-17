%module pbconsensus

%{
#define SWIG_FILE_WITH_INIT
%}

%pythoncode
%{
__version__ = "0.9.0"
%}

%ignore *::operator[];
%ignore *::operator=;
%ignore *::operator<<;
%ignore *::operator std::string;

%ignore boost::noncopyable;
namespace boost { class noncopyable {}; }

%include <std_string.i>
%include <std_vector.i>
%include <stdint.i>

#ifdef SWIGPYTHON
%include "numpy.i"
#endif

namespace std
{
    %template(Uint8Vector) std::vector<uint8_t>;
}

// no include dependencies
%include "Exceptions.i"
%include "ModelConfig.i"
%include "Mutation.i"
%include "Polish.i"
%include "Read.i"

// after Read.i
%include "Integrator.i"
%include "Poa.i"

%init
%{
#ifdef SWIGPYTHON
  import_array();
#endif
%}
