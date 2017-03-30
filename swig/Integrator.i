
%{
#include <pacbio/consensus/AbstractIntegrator.h>
#include <pacbio/consensus/Integrator.h>
%}

%feature("notabstract") PacBio::Consensus::Integrator;

py_tp_str(PacBio::Consensus::AbstractIntegrator);
py_tp_str(PacBio::Consensus::Integrator);

%include <pacbio/consensus/AbstractIntegrator.h>
%include <pacbio/consensus/Integrator.h>
