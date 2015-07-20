//
//  TransitionParameters.cpp
//  ConsensusCoreXCode
//
//  Created by Nigel Delaney on 2/23/15.
//  Copyright (c) 2015 Pacific Biosciences. All rights reserved.
//

#include <ConsensusCore/Arrow/TransitionParameters.hpp>
#include <ConsensusCore/Arrow/MathUtils.hpp>

namespace ConsensusCore {
namespace Arrow {

    TransitionParameters::TransitionParameters(double match, double stick,
                                               double branch, double deletion) :
                                               Match(match), Stick(stick),
                                               Branch(branch),Deletion(deletion) {}


    TransitionParameters::TransitionParameters() :
    Match(0.0) ,Stick(0.0), Branch(0.0), Deletion(0.0) {}

    TransitionParameters::TransitionParameters(const TransitionParameters &other) :
    Match(other.Match), Stick(other.Stick), Branch(other.Branch), Deletion(other.Deletion) {}

    double TransitionParameters::CalculateTotal() const {
        return Match + Stick + Branch + Deletion; //logsumlog(Match, Stick, Branch, Deletion);

    }

    void TransitionParameters::RemoveConstant(double value) {
        Match /= value;
        Stick /= value;
        Branch /= value;
        Deletion /= value;
    }
}
}
