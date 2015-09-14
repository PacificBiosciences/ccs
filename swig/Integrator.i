
%{
#include <pacbio/consensus/Integrator.h>
using namespace PacBio::Consensus;
%}

%feature("notabstract") MonoMolecularIntegrator;
%feature("notabstract") MultiMolecularIntegrator;

%include <pacbio/consensus/Integrator.h>
