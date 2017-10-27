// TODO(lhepler): remove these, figure out %ignore
%warnfilter(509) PacBio::Consensus::Mutation;

%{
#include <pacbio/consensus/Mutation.h>
%}

py_tp_str(PacBio::Consensus::Mutation);
py_tp_str(PacBio::Consensus::ScoredMutation);

%include <pacbio/consensus/Mutation.h>

namespace std {
    %ignore vector<PacBio::Consensus::Mutation>::vector(size_type);
    %ignore vector<PacBio::Consensus::Mutation>::resize;
    %template(MutationVector) vector<PacBio::Consensus::Mutation>;

    %ignore vector<PacBio::Consensus::ScoredMutation>::vector(size_type);
    %ignore vector<PacBio::Consensus::ScoredMutation>::resize;
    %template(ScoredMutationVector) vector<PacBio::Consensus::ScoredMutation>;
 }
