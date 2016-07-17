
%{
#include <pacbio/align/AlignConfig.h>
#include <pacbio/align/PairwiseAlignment.h>
#include <pacbio/align/AffineAlignment.h>
#include <pacbio/align/LinearAlignment.h>
%}

%newobject PacBio::Align::Align;
%newobject PacBio::Align::AlignAffine;
%newobject PacBio::Align::AlignAffineIupac;
%newobject PacBio::Align::AlignLinear;

%include <pacbio/align/AlignConfig.h>
%include <pacbio/align/PairwiseAlignment.h>
%include <pacbio/align/AffineAlignment.h>
%include <pacbio/align/LinearAlignment.h>
