//
//  ContextParameterProvider.cpp
//  ConsensusCoreXCode
//
//  Created by Nigel Delaney on 2/27/15.
//  Copyright (c) 2015 Pacific Biosciences. All rights reserved.
//

#include <ConsensusCore/Arrow/ContextParameterProvider.hpp>

using namespace std;

namespace ConsensusCore {
namespace Arrow {
    
    SNR::SNR(double a, double c, double g, double t)
        : A{a}, C{c}, G{g}, T{t}
    { }

    //Rows are Dark, Match, Stick (Branch is the reference)
    //Columns are Intercept, SNR, SNR^2, SNR^3
    // Fit for context:  AA
    Matrix<double>  AA  = {
        { 3.76122480667588, -0.536010820176981, 0.0275375059387171, -0.000470200724345621  },
        { 3.57517725358548, -0.0257545295375707, -0.000163673803286944, 5.3256984681724e-06  },
        { 0.858421613302247, -0.0276654216841666, -8.85549766507732e-05, -4.85355908595337e-05  } };
    // Fit for context:  CC
    Matrix<double>  CC  = {
        { 5.66725538674764, -1.10462196933913, 0.0879811093908922, -0.00259393800835979  },
        { 4.11682756767018, -0.124758322644639, 0.00659795177909886, -0.000361914629195461  },
        { 3.17103818507405, -0.729020290806687, 0.0749784690396837, -0.00262779517495421  } };
    // Fit for context:  GG
    Matrix<double>  GG  = {
        { 3.81920778703052, -0.540309003502589, 0.0389569264893982, -0.000901245733796236  },
        { 3.31322216145728, 0.123514009118836, -0.00807401406655071, 0.000230843924466035  },
        { 2.06006877520527, -0.451486652688621, 0.0375212898173045, -0.000937676250926241  } };
    // Fit for context:  TT
    Matrix<double>  TT  = {
        { 5.39308368236762, -1.32931568057267, 0.107844580241936, -0.00316462903462847  },
        { 4.21031404956015, -0.347546363361823, 0.0293839179303896, -0.000893802212450644  },
        { 2.33143889851302, -0.586068444099136, 0.040044954697795, -0.000957298861394191  } };
    // Fit for context:  NA
    Matrix<double>  NA  = {
        { 2.35936060895653, -0.463630601682986, 0.0179206897766131, -0.000230839937063052  },
        { 3.22847830625841, -0.0886820214931539, 0.00555981712798726, -0.000137686231186054  },
        { -0.101031042923432, -0.0138783767832632, -0.00153408019582419, 7.66780338484727e-06  } };
    // Fit for context:  NC
    Matrix<double>  NC  = {
        { 5.956054206161, -1.71886470811695, 0.153315470604752, -0.00474488595513198  },
        { 3.89418464416296, -0.174182841558867, 0.0171719290275442, -0.000653629721359769  },
        { 2.40532887070852, -0.652606650098156, 0.0688783864119339, -0.00246479494650594  } };
    // Fit for context:  NG
    Matrix<double>  NG  = {
        { 3.53508304630569, -0.788027301381263, 0.0469367803413207, -0.00106221924705805  },
        { 2.85440184222226, 0.166346531056167, -0.0166161828155307, 0.000439492705370092  },
        { 0.238188180807376, 0.0589443522886522, -0.0123401045958974, 0.000336854126836293  } };
    // Fit for context:  NT
    Matrix<double>  NT  = {
        { 5.36199280681367, -1.46099908985536, 0.126755291030074, -0.0039102734460725  },
        { 3.41597143103046, -0.066984162951578, 0.0138944877787003, -0.000558939998921912  },
        { 1.37371376794871, -0.246963827944892, 0.0209674231346363, -0.000684856715039738  } };
    
    
    static std::unordered_map<std::string, Matrix<double>* > parameter_store = { {"AA", &AA},{"CC", &CC},{"GG", &GG},{"NA", &NA},{"NC", &NC},{"NG", &NG},{"NT", &NT},{"TT", &TT}};
    
    TransitionParameters
    ContextParameterProvider::GetTransitionParameters(const string& context, const SNR& snrs)
    {
        auto params = *parameter_store[context];
        //Get the snr for the relevant channel
        auto channel  = context.at(1);
        double snr;
        switch(channel) {
                case 'A':
                    snr = snrs.A;
                    break;
                case 'C':
                    snr = snrs.C;
                    break;
                case 'G':
                    snr = snrs.G;
                    break;
                case 'T':
                    snr = snrs.T;
                    break;
            default:
                throw;
        }
        double snr2 = snr * snr;
        double snr3 = snr2 * snr;
        
        double predicts[3]; // Represents the XB portion
        double sum = 1.0;
        // Calculate each values contribution
        for(int i=0; i< 3; i ++) {
            auto xb = params[i][0] + snr * params[i][1] + snr2 * params[i][2] + snr3 * params[i][3];
            xb = exp(xb);
            predicts[i] = xb;
            sum += xb;
        }
        
        double branch = 1.0 / sum; // match probability is the reference, or 1 / sum
        
        // Now get the probabilities
        for(int i=0; i< 3; i++) {
            predicts[i] = predicts[i] / sum;
        }
        TransitionParameters tp(predicts[1], predicts[2], branch, predicts[0]);        
        return tp;
    }
    
}
}
