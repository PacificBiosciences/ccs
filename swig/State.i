%{
#include <pacbio/data/State.h>
%}

%ignore PacBio::Data::StateName;
%include <pacbio/data/State.h>

namespace std {
    %template(StateVector) vector<PacBio::Data::State>;
}
