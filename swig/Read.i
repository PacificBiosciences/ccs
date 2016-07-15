
// TODO(lhepler): remove these, figure out %ignore
%warnfilter(509) PacBio::Consensus::Read::Read;
%warnfilter(509) PacBio::Consensus::MappedRead::MappedRead;

%{
#include <pacbio/consensus/StrandType.h>
#include <pacbio/consensus/Read.h>
%}

%feature("notabstract") PacBio::Consensus::Read;
%feature("notabstract") PacBio::Consensus::MappedRead;

%include <pacbio/consensus/StrandType.h>
%include <pacbio/consensus/Read.h>
