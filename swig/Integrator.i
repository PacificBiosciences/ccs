
%{
#include <pacbio/consensus/Integrator.h>
%}

%feature("notabstract") PacBio::Consensus::Integrator;

py_tp_str(PacBio::Consensus::Integrator);

%include <pacbio/consensus/Integrator.h>
