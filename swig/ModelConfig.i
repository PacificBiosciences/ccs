
%{
#include <pacbio/consensus/ModelConfig.h>
using namespace PacBio::Consensus;
%}

%ignore PacBio::Consensus::TemplatePosition;
%ignore PacBio::Consensus::MoveType;
%ignore PacBio::Consensus::ModelConfig;

%include <pacbio/consensus/ModelConfig.h>
