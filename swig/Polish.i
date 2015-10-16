
%{
#include <pacbio/consensus/Polish.h>
%}

// for results from Polish
%typemap(out) std::tuple<bool, size_t, size_t> {
    $result = PyTuple_Pack(3, PyBool_FromLong(std::get<0>($1)),
                              PyInt_FromSize_t(std::get<1>($1)),
                              PyInt_FromSize_t(std::get<2>($1)));
}

%include <pacbio/consensus/Polish.h>
