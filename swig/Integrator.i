
%{
#include <pacbio/consensus/Integrator.h>
%}

%ignore PacBio::Consensus::AbstractIntegrator::ReadState;

%feature("notabstract") PacBio::Consensus::MonoMolecularIntegrator;
%feature("notabstract") PacBio::Consensus::MultiMolecularIntegrator;

py_tp_str(PacBio::Consensus::AbstractIntegrator);
py_tp_str(PacBio::Consensus::MonoMolecularIntegrator);
py_tp_str(PacBio::Consensus::MultiMolecularIntegrator);

%include <pacbio/consensus/Integrator.h>
