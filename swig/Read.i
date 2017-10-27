
// TODO(lhepler): remove these, figure out %ignore
//%warnfilter(503) PacBio::Data::Read;
%warnfilter(509) PacBio::Data::Read;
%warnfilter(509) PacBio::Data::MappedRead;

%{
#include <pacbio/data/StrandType.h>
#include <pacbio/data/Read.h>
%}

%feature("notabstract") PacBio::Data::Read;
%feature("notabstract") PacBio::Data::MappedRead;

%ignore operator std::vector<float>;

%include <pacbio/data/StrandType.h>
%include <pacbio/data/Read.h>
