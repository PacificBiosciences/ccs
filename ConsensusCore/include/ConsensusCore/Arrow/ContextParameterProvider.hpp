//
//  ContextParameterProvider.h
//  ConsensusCoreXCode
//
//  Created by Nigel Delaney on 2/27/15.
//  Copyright (c) 2015 Pacific Biosciences. All rights reserved.
//

#pragma once

#include <stdio.h>
#include <unordered_map>
#include <boost/numeric/ublas/matrix.hpp>
#include <cmath>
#include <vector>

#include <ConsensusCore/Arrow/TransitionParameters.hpp>

namespace ConsensusCore {
namespace Arrow {

    template<class T>
    using Matrix = std::vector<std::vector<T> >;

    /**
     A data structure to store SNR values.
     */
    class SNR {
    public:
        const double A, C, G, T;
        SNR(double A, double C, double G, double T);

        // support arbitrary loading from vector
        template<typename RealType>
        SNR(const std::vector<RealType>& snrs)
            : A{snrs[0]}
            , C{snrs[1]}
            , G{snrs[2]}
            , T{snrs[3]}
        { }
    };
    
    /**
     
     This class is designed to provide the relative transition probabilities for
     a given dinculeotide context at a given SNR value.
     @param <#parameter#>
     @returns <#retval#>
     @exception <#throws#>
     */
    class ContextParameterProvider {
        public:
        static TransitionParameters GetTransitionParameters(const std::string& context, const SNR& snrs);
    };
}
}
