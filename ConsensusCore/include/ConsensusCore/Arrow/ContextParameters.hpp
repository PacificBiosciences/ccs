//
//  ContextParameters.h
//  ConsensusCoreXCode
//
//  Created by Nigel Delaney on 3/1/15.
//  Copyright (c) 2015 Pacific Biosciences. All rights reserved.
//

#pragma once

#include <stdio.h>
#include <unordered_map>
#include <string>
#include <vector>

#include <ConsensusCore/Arrow/ContextParameterProvider.hpp>

namespace ConsensusCore {
namespace Arrow {
    /**
     This class represents a collection of context parameters for a given set of
     SNR values.  It stores the transition parameters for the 
     @param <#parameter#>
     @returns <#retval#>
     @exception <#throws#>
     */
    class ContextParameters {
    public:
        ContextParameters(SNR snrs);
        // Empty constructor for swig
        ContextParameters();
        // Copy constructor for swig
        ContextParameters(const ContextParameters& arg);
        TransitionParameters GetParametersForContext(char bp1, char bp2) const;
    private:
        std::unordered_map<std::string, TransitionParameters> param_map;
        std::vector<std::string> contexts {"AA", "NA", "CC", "NC", "TT", "NT", "GG", "NG"};
    };
}
}
