
%{
#include <pacbio/consensus/ModelConfig.h>
#include <pacbio/consensus/ModelSelection.h>
%}

%ignore PacBio::Consensus::TemplatePosition;
%ignore PacBio::Consensus::MoveType;
%ignore PacBio::Consensus::ModelConfig;

%include <pacbio/consensus/ModelConfig.h>
%include <pacbio/consensus/ModelSelection.h>
