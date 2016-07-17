
// TODO(lhepler): remove these, figure out %ignore
%warnfilter(509) PacBio::Data::Read::Read;
%warnfilter(509) PacBio::Data::MappedRead::MappedRead;

%{
#include <pacbio/data/StrandType.h>
#include <pacbio/data/Read.h>
%}

%feature("notabstract") PacBio::Data::Read;
%feature("notabstract") PacBio::Data::MappedRead;

%include <pacbio/data/StrandType.h>
%include <pacbio/data/Read.h>
