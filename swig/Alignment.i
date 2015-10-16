
%{
#include <pacbio/consensus/align/AlignConfig.h>
#include <pacbio/consensus/align/PairwiseAlignment.h>
#include <pacbio/consensus/align/AffineAlignment.h>
#include <pacbio/consensus/align/LinearAlignment.h>
%}

%newobject PacBio::Consensus::Align;
%newobject PacBio::Consensus::AlignAffine;
%newobject PacBio::Consensus::AlignAffineIupac;
%newobject PacBio::Consensus::AlignLinear;

%include <pacbio/consensus/align/AlignConfig.h>
%include <pacbio/consensus/align/PairwiseAlignment.h>
%include <pacbio/consensus/align/AffineAlignment.h>
%include <pacbio/consensus/align/LinearAlignment.h>
