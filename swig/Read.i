
%{
#include <pacbio/consensus/Read.h>
using namespace PacBio::Consensus;
%}

%feature("notabstract") Read;
%feature("notabstract") MappedRead;

%include <pacbio/consensus/Read.h>
