
%{
#include <pacbio/consensus/Template.h>
%}

%feature("notabstract") PacBio::Consensus::Template;

%ignore PacBio::Consensus::AbstractTemplate;
%ignore PacBio::Consensus::VirtualTemplate;

%include <pacbio/consensus/Template.h>
