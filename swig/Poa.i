
%{
#include <pacbio/denovo/PoaGraph.h>
#include <pacbio/denovo/PoaConsensus.h>
%}

%include <pacbio/denovo/PoaGraph.h>

%newobject PacBio::Poa::PoaConsensus::FindConsensus;

%include <pacbio/denovo/PoaConsensus.h>
