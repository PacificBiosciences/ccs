
%{
#include <pacbio/consensus/Read.h>
%}

%feature("notabstract") PacBio::Consensus::Read;
%feature("notabstract") PacBio::Consensus::MappedRead;

%include <pacbio/consensus/Read.h>
