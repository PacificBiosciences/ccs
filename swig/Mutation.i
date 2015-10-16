
%{
#include <pacbio/consensus/Mutation.h>
%}

py_tp_str(PacBio::Consensus::Mutation);
py_tp_str(PacBio::Consensus::ScoredMutation);

%include <pacbio/consensus/Mutation.h>
