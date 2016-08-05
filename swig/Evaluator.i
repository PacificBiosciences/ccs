%{
#include <pacbio/consensus/Evaluator.h>
#include <pacbio/consensus/EasyReadScorer.h>
%}

%newobject PacBio::Consensus::EasyReadScorer::MakeEvaluator;

%include <pacbio/consensus/Evaluator.h>
%include <pacbio/consensus/EasyReadScorer.h>
