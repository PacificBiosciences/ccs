%module pbconsensus

%{
#define SWIG_FILE_WITH_INIT
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

%include "Exceptions.i"
%include "Integrator.i"
%include "Mutation.i"
%include "Poa.i"
%include "Polish.i"
%include "Read.i"

%init
%{
#ifdef SWIGPYTHON
  import_array();
#endif
%}
