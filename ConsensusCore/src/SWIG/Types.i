//
// Include this file *first* to pull in forward declarations.
//

%{
/* Includes the header in the wrapper code */
#include <ConsensusCore/Interval.hpp>
#include <ConsensusCore/Types.hpp>
using namespace ConsensusCore;
%}

%include <std_pair.i>
%include <std_string.i>
%include <std_vector.i>
%include <std_map.i>
%include "stdint.i"

%include <ConsensusCore/Types.hpp>

// Need to see this definition before we see intervals.
namespace std {
  %template(IntPair)                std::pair<int, int>;
};

%include <ConsensusCore/Interval.hpp>

namespace std {
  %template(IntervalVector)         std::vector<ConsensusCore::Interval>;
  %template(IntVector)              std::vector<int>;
  %template(FloatVector)            std::vector<float>;
  %template(StringVector)           std::vector<string>;
  %template(FeaturesVector)         std::vector<const ConsensusCore::SequenceFeatures*>;
};
