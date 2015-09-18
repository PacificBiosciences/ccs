
%{
#include <pacbio/consensus/Mutation.h>
using namespace PacBio::Consensus;
%}

py_tp_str(PacBio::Consensus::Mutation);
py_tp_str(PacBio::Consensus::ScoredMutation);

%include <pacbio/consensus/Mutation.h>
