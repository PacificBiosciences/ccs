
%{
#include <pacbio/consensus/poa/AlignConfig.h>
#include <pacbio/consensus/poa/PoaGraph.h>
#include <pacbio/consensus/poa/PoaConsensus.h>
using namespace PacBio::Consensus;
%}

%include <pacbio/consensus/poa/AlignConfig.h>
%include <pacbio/consensus/poa/PoaGraph.h>

%newobject PacBio::Consensus::PoaConsensus::FindConsensus;

%include <pacbio/consensus/poa/PoaConsensus.h>
