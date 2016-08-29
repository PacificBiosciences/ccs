
%{
#include <pacbio/data/State.h>
%}

%include <pacbio/data/State.h>

namespace std {
    %template(StateVector) vector<PacBio::Data::State>;
}
