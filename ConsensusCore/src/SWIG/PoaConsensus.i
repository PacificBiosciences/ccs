
%{
/* Includes the header in the wrapper code */
#include <ConsensusCore/Poa/PoaGraph.hpp>
#include <ConsensusCore/Poa/PoaConsensus.hpp>
using namespace ConsensusCore;
%}


#ifdef SWIGCSHARP
 // SWIG needs this to properly override the ToString method
%csmethodmodifiers ConsensusCore::PoaConsensus::ToString() const "public override"
#endif // SWIGCSHARP

%include <ConsensusCore/Poa/PoaGraph.hpp>

%newobject ConsensusCore::PoaConsensus::FindConsensus;

%include <ConsensusCore/Poa/PoaConsensus.hpp>
