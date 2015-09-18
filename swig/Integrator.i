
%{
#include <pacbio/consensus/Integrator.h>
using namespace PacBio::Consensus;
%}

%feature("notabstract") MonoMolecularIntegrator;
%feature("notabstract") MultiMolecularIntegrator;

py_tp_str(PacBio::Consensus::AbstractIntegrator);
py_tp_str(PacBio::Consensus::MonoMolecularIntegrator);
py_tp_str(PacBio::Consensus::MultiMolecularIntegrator);

%include <pacbio/consensus/Integrator.h>
