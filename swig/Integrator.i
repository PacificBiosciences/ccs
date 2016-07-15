
%{
#include <pacbio/consensus/State.h>
#include <pacbio/consensus/AbstractIntegrator.h>
#include <pacbio/consensus/MonoMolecularIntegrator.h>
#include <pacbio/consensus/MultiMolecularIntegrator.h>
%}

%feature("notabstract") PacBio::Consensus::MonoMolecularIntegrator;
%feature("notabstract") PacBio::Consensus::MultiMolecularIntegrator;

py_tp_str(PacBio::Consensus::AbstractIntegrator);
py_tp_str(PacBio::Consensus::MonoMolecularIntegrator);
py_tp_str(PacBio::Consensus::MultiMolecularIntegrator);

%include <pacbio/consensus/State.h>
%include <pacbio/consensus/AbstractIntegrator.h>
%include <pacbio/consensus/MonoMolecularIntegrator.h>
%include <pacbio/consensus/MultiMolecularIntegrator.h>
