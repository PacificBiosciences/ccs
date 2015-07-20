%{
#include <ConsensusCore/Mutation.hpp>
#include <ConsensusCore/Quiver/MutationEnumerator.hpp>
using namespace ConsensusCore;
%}

// SWIG needs this to properly override the ToString method
#ifdef SWIGCSHARP
%csmethodmodifiers ConsensusCore::Mutation::ToString() const "public override"
#endif // SWIGCSHARP

 // SWIG now seems to be incorrectly deciding that these are
 // abstract classes, so we have to tell it otherwise
%feature("notabstract") AllSingleBaseMutationEnumerator;
%feature("notabstract") UniqueSingleBaseMutationEnumerator;
%feature("notabstract") RepeatMutationEnumerator;
%feature("notabstract") DinucleotideRepeatMutationEnumerator;

%newobject *::Mutations();
%newobject *::Mutations(int, int);

%include <std_pair.i>
%include <std_vector.i>

namespace std {
    %template(MutationVector)        std::vector<ConsensusCore::Mutation>;
    %template(ScoredMutationVector)  std::vector<ConsensusCore::ScoredMutation>;
};

%include <ConsensusCore/Mutation.hpp>
%include <ConsensusCore/Quiver/MutationEnumerator.hpp>
