%{
/* Includes the header in the wrapper code */
#include <ConsensusCore/Align/AlignConfig.hpp>
#include <ConsensusCore/Align/PairwiseAlignment.hpp>
#include <ConsensusCore/Align/AffineAlignment.hpp>
#include <ConsensusCore/Align/LinearAlignment.hpp>
using namespace ConsensusCore;
%}

%include "Types.i"

%newobject Align;
%newobject AlignAffine;
%newobject AlignAffineIupac;
%newobject AlignLinear;

%include <ConsensusCore/Align/AlignConfig.hpp>
%include <ConsensusCore/Align/PairwiseAlignment.hpp>
%include <ConsensusCore/Align/AffineAlignment.hpp>
%include <ConsensusCore/Align/LinearAlignment.hpp>
