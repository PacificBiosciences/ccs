
%{
#include <pacbio/UnanimityVersion.h>
%}

%pythoncode
%{
__version__ = UnanimityVersion()
%}

%include <pacbio/UnanimityVersion.h>
