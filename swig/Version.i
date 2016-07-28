
%{
#include <pacbio/Version.h>
%}

%pythoncode
%{
__version__ = UnanimityVersion()
%}

%include <pacbio/Version.h>
