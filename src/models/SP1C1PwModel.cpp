// Copyright (c) 2011-2016, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

#include <cassert>
#include <cmath>
#include <memory>
#include <stdexcept>

#include <pacbio/consensus/ModelConfig.h>
#include <pacbio/consensus/Read.h>

#include "../ModelFactory.h"
#include "../Recursor.h"

namespace PacBio {
namespace Consensus {
namespace {

constexpr double kCounterWeight = 15.0;
const size_t OUTCOME_NUMBER = 12;
const size_t CONTEXT_NUMBER = 16;
class SP1C1PwModel : public ModelConfig
{
    REGISTER_MODEL(SP1C1PwModel);

public:
    static std::set<std::string> Names() { return {"S/P1-C1.1"}; }
    SP1C1PwModel(const SNR& snr);
    std::unique_ptr<AbstractRecursor> CreateRecursor(std::unique_ptr<AbstractTemplate>&& tpl,
                                                     const MappedRead& mr, double scoreDiff) const;
    std::vector<TemplatePosition> Populate(const std::string& tpl) const;
    double ExpectedLogLikelihoodForMatchEmission(uint8_t prev, uint8_t curr, bool secondMoment) const;
    double ExpectedLogLikelihoodForStickEmission(uint8_t prev, uint8_t curr, bool secondMoment) const;
    double ExpectedLogLikelihoodForBranchEmission(uint8_t prev, uint8_t curr, bool secondMoment) const;
    


private:
    SNR snr_;
    double ctxTrans_[CONTEXT_NUMBER][4];
    double cachedEmissionExpectations_[CONTEXT_NUMBER][3][2];
    double ExpectedLogLikelihoodOfOutcomeRow(const int index, const uint8_t prev, const uint8_t curr, const bool secondMoment) const;
};

REGISTER_MODEL_IMPL(SP1C1PwModel);

// TODO(lhepler) comments regarding the CRTP
class SP1C1PwRecursor : public Recursor<SP1C1PwRecursor>
{
public:
    SP1C1PwRecursor(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr,
                    double scoreDiff);
    static inline std::vector<uint8_t> EncodeRead(const MappedRead& read);
    static inline double EmissionPr(MoveType move, uint8_t emission, uint8_t prev, uint8_t curr);
    virtual double UndoCounterWeights(size_t nEmissions) const;
};

    
constexpr double emissionPmf[3][CONTEXT_NUMBER][OUTCOME_NUMBER] = {
    {// matchPmf
        {   0.0516695237,  0.000856759021,  0.000805765869,  6.40098033e-05,    0.0735681929,  0.000490647727,   0.00171763848,  1.66588276e-05,     0.864949369,   0.00298262543,    0.0027926093,  6.28389773e-05},
        {  0.00543436007,     0.011889156,  0.000101241451,  0.000115879948,    0.0161576112,    0.0379073447,   0.00013138862,  0.000162963331,    0.0238359857,     0.903407615,  0.000310425863,   0.00051517498},
        {  0.00078327569,   8.9195486e-05,    0.0213823858,   0.00233242884,   0.00122583139,  2.26787869e-05,    0.0797753558,    0.0023406575,  0.000536349203,  0.000218984243,     0.881882162,   0.00937678742},
        { 0.000153647232,   0.00010627987,    0.0103924584,    0.0168897121,  0.000528420929,  0.000114998697,     0.034445993,    0.0461495153,  0.000202642175,  0.000688915405,    0.0797339416,     0.810568344},
        {   0.0217684866,  0.000615574217,   0.00123796174,  0.000194115893,    0.0630331228,  0.000312092976,   0.00217859654,  1.70111886e-05,     0.905859233,   0.00150065174,    0.0030195303,   0.00023882051},
        {    0.043036651,    0.0109154406,  0.000137917521,  6.63118913e-05,    0.0135411051,    0.0336696383,  3.68321175e-05,   9.7192835e-05,    0.0198015123,     0.878552904,  1.95554178e-05,  9.54805479e-05},
        {  0.00018454615,  1.51427373e-05,   0.00855522704,   0.00108940771,  0.000328562866,  1.41858457e-05,    0.0341491789,   0.00128003844,   0.00124774351,  0.000256778116,     0.944704089,   0.00815078658},
        { 0.000244636834,  5.95364998e-05,    0.0103025634,    0.0176534311,  0.000417091732,  3.89027835e-05,    0.0347137702,    0.0495173978,  0.000658139356,  0.000481549284,    0.0798548882,      0.80602469},
        {   0.0120953705,  0.000297201253,  0.000336831573,  0.000109810239,    0.0339973245,  0.000295424749,  0.000821327226,   7.0554379e-05,     0.947447668,   0.00298212225,   0.00115060424,  0.000366663193},
        {  0.00202072785,    0.0061102332,  0.000138730316,  8.82424379e-05,   0.00783346136,    0.0193899073,  0.000128803697,   5.1853045e-05,    0.0145084431,     0.948899725,  0.000308586587,   0.00049984933},
        {   0.0362076945,  3.82567383e-05,   0.00699047085,  0.000843806016,   0.00025900091,  2.74289133e-05,    0.0275479876,   0.00133679409,   0.00089127495,  0.000135756189,     0.918742407,   0.00694931573},
        { 0.000148143946,  5.09628138e-05,   0.00492091587,   0.00825067123,  3.34491844e-05,  1.52215475e-05,    0.0185009947,    0.0282098475,     0.000211194,  0.000843725462,    0.0534485418,     0.885335084},
        {   0.0233179774,  0.000994145372,  0.000594429081,  3.30759279e-05,     0.061000414,   0.00042148067,   0.00179973978,     7.98899e-05,     0.906095769,   0.00277286156,   0.00269906625,  0.000153458766},
        {  0.00252578813,    0.0058253842,  2.66040065e-05,  8.55149972e-05,   0.00783559799,     0.020356316,  0.000125999422,  9.67021798e-05,    0.0145272461,     0.948105104,  0.000231890377,  0.000228741056},
        { 0.000344291837,  1.71744233e-05,    0.0093318888,   0.00101109441,  0.000518176918,  2.49752238e-05,    0.0364757971,   0.00118973194,   0.00131474596,  0.000128603181,     0.942187364,   0.00743159244},
        {   0.0286354873,  0.000348696974,   0.00669477573,   0.00817040337,   0.00263533574,  0.000832163375,    0.0187454532,    0.0241172554,    0.0272027554,    0.0261244204,    0.0786962737,     0.777776291}},

    {// branchPmf
        {    0.337545783,  0.000189833088,  0.000189833088,  0.000189833088,     0.167897832,  0.000189833088,  0.000189833088,  0.000189833088,     0.491898721,  0.000189833088,  0.000189833088,  0.000189833088},
        { 0.000455211501,      0.10721842,  0.000455211501,  0.000455211501,  0.000455211501,    0.0880290507,  0.000455211501,  0.000455211501,  0.000455211501,     0.798379568,  0.000455211501,  0.000455211501},
        { 0.000139076929,  0.000139076929,     0.277526154,  0.000139076929,  0.000139076929,  0.000139076929,     0.164738707,  0.000139076929,  0.000139076929,  0.000139076929,     0.555788062,  0.000139076929},
        { 0.000340988677,  0.000340988677,  0.000340988677,     0.107667327,  0.000340988677,  0.000340988677,  0.000340988677,      0.11481587,  0.000340988677,  0.000340988677,  0.000340988677,     0.772742962},
        {     0.15792522,  4.71456125e-05,  4.71456125e-05,  4.71456125e-05,     0.141752834,  4.71456125e-05,  4.71456125e-05,  4.71456125e-05,     0.699661908,  4.71456125e-05,  4.71456125e-05,  4.71456125e-05},
        { 8.48749074e-05,    0.0356794691,  8.47792483e-05,  8.47792483e-05,  8.50080149e-05,    0.0463288295,  8.47792483e-05,  8.47792483e-05,  8.85882375e-05,     0.916800659,  8.47792483e-05,  8.47792483e-05},
        { 7.09844065e-05,  7.09844065e-05,     0.228163146,  7.09844065e-05,  7.09844065e-05,  7.09844065e-05,     0.125258165,  7.09844065e-05,  7.09844065e-05,  7.09844065e-05,     0.645584908,  7.09844065e-05},
        { 0.000294417195,  0.000294417195,  0.000294417195,    0.0762864785,  0.000294417195,  0.000294417195,  0.000294417195,     0.115945928,  0.000294417195,  0.000294417195,  0.000294417195,     0.803645753},
        {    0.116533833,  5.79322675e-05,  5.79322675e-05,  5.79322675e-05,    0.0942835568,  5.79322675e-05,  5.79322675e-05,  5.79322675e-05,     0.788371558,  5.79322675e-05,  5.79322675e-05,  5.79322675e-05},
        { 0.000156684083,    0.0527691682,  0.000156684083,  0.000156684083,  0.000156684083,    0.0520229677,  0.000156684083,  0.000156684083,  0.000156684083,     0.893014287,  0.000156684083,  0.000156684083},
        { 0.000167018487,  0.000165322717,     0.312421562,  0.000165322717,  0.000165401501,  0.000165322717,     0.133726114,  0.000165322717,  0.000165734487,  0.000165322717,      0.55153562,  0.000165322717},
        { 0.000203462057,  0.000203462057,  0.000203462057,    0.0471982179,  0.000203462057,  0.000203462057,  0.000203462057,     0.049201974,  0.000203462057,  0.000203462057,  0.000203462057,     0.900751339},
        {    0.119585426,  6.14072028e-05,  6.14072028e-05,  6.14072028e-05,     0.116129123,  6.14072028e-05,  6.14072028e-05,  6.14072028e-05,      0.76342575,  6.14072028e-05,  6.14072028e-05,  6.14072028e-05},
        { 0.000110720498,    0.0349410744,  0.000110720498,  0.000110720498,  0.000110720498,    0.0453942421,  0.000110720498,  0.000110720498,  0.000110720498,     0.918114596,  0.000110720498,  0.000110720498},
        { 7.09455464e-05,  7.09455464e-05,     0.186570455,  7.09455464e-05,  7.09455464e-05,  7.09455464e-05,      0.11383208,  7.09455464e-05,  7.09455464e-05,  7.09455464e-05,     0.698604227,  7.09455464e-05},
        { 5.44400066e-05,  4.58781747e-05,  4.58781747e-05,    0.0513058841,  4.58791482e-05,  4.58781747e-05,  4.58781747e-05,     0.081923538,  4.95222541e-05,  4.58781747e-05,  4.58781747e-05,     0.866116077}},

    {// stickPmf
        { 0.000154038684,    0.0291154091,     0.423444453,    0.0352058932,  0.000154038684,    0.0176086145,     0.149856392,    0.0151919119,  0.000154038684,    0.0916984378,     0.168808535,     0.067838043},
        {     0.13839578,  9.92440537e-05,     0.329722566,    0.0318555909,    0.0371088704,  9.92440537e-05,     0.148457457,    0.0195393534,     0.112987247,  9.92440537e-05,     0.103583418,    0.0775557653},
        {     0.39242986,    0.0361688498,   0.00033884198,     0.138030777,    0.0355767848,    0.0102554086,   0.00033884198,     0.054319306,    0.0665562733,    0.0785778949,   0.00033884198,     0.185374109},
        {    0.265583401,    0.0323250154,     0.267620652,  0.000144990997,    0.0553401869,    0.0100214732,      0.12586624,  0.000144990997,     0.081077577,    0.0457693156,     0.115236212,  0.000144990997},
        { 8.86069101e-05,    0.0224045604,     0.391863179,    0.0374698777,  8.86069101e-05,     0.014879736,     0.127064528,    0.0148023081,  8.86069101e-05,     0.199043674,     0.127110897,    0.0646523839},
        {    0.132279199,  5.64661147e-05,     0.321806183,    0.0323685244,    0.0725530116,  4.63524839e-05,     0.160579431,    0.0137370025,     0.108005126,  4.69523474e-05,     0.116162938,    0.0421336554},
        {    0.320519663,    0.0239698473,   0.00011228542,    0.0412339441,     0.124873594,    0.0131178903,   0.00011228542,    0.0214946836,     0.139243924,     0.223216618,   0.00011228542,     0.091431553},
        {    0.236440465,    0.0125729978,     0.174982446,  0.000114967022,     0.072742801,    0.0157082833,     0.107504931,  0.000114967022,    0.0812575292,     0.199021401,    0.0988494102,  0.000114967022},
        { 0.000172377278,    0.0252780963,     0.441347929,    0.0431573519,  0.000172377278,    0.0109779693,     0.122083706,      0.01370133,  0.000172377278,     0.112156842,     0.126224513,     0.103693244},
        {     0.11616354,   4.4550521e-05,      0.38730196,    0.0316267288,     0.048597637,   4.4550521e-05,     0.217316542,    0.0147582454,    0.0443979757,   4.4550521e-05,    0.0941246042,    0.0453563631},
        {    0.351437251,    0.0237786623,   0.00123104756,    0.0609100017,     0.132002208,    0.0127583294,  0.000288446433,    0.0203701885,       0.1316581,     0.080542757,  0.000252892079,      0.18353973},
        {    0.298492049,    0.0269851185,      0.29827165,  0.000200099143,     0.102490322,    0.0117983684,     0.121179733,  0.000200099143,    0.0770158895,    0.0380659958,    0.0241000812,  0.000200099143},
        { 0.000141627536,      0.01777018,     0.301881108,    0.0226933622,  0.000141627536,    0.0109104656,     0.114187821,    0.0163453226,  0.000141627536,     0.102595478,     0.208037711,     0.204445531},
        {    0.132859004,  7.16119861e-05,       0.2406106,    0.0231724181,     0.057038948,  7.16119861e-05,    0.0925183577,    0.0141944935,    0.0659598334,  7.16119861e-05,     0.125422694,     0.247650755},
        {    0.355614563,    0.0234102543,  0.000187722927,    0.0564654303,     0.109173536,    0.0154941999,  0.000187722927,    0.0155496293,     0.106690121,    0.0539175558,  0.000187722927,     0.262182928},
        {    0.152511749,    0.0105044023,     0.228700197,  9.13331245e-05,    0.0553648837,   0.00620728656,     0.191510525,  7.33906221e-05,     0.048230451,    0.0444486569,     0.261928846,  7.37812688e-05}}};

constexpr double transProbs[16][3][4] = {
    // Fit for context:  AA
    {
        { -11.8402163277545, 4.12000773672109, -0.661725654820908, 0.0337357107347639  },
        { -3.10509518253287, 0.275600076181856, -0.101980967380149, 0.00728037713082447  },
        { 8.79817734479091, -4.80619827278369, 0.668699377698361, -0.0323014351665758  }
    },
    // Fit for context:  AC
    {
        { 3.76105501880563, -3.52354687054574, 0.518155785867385, -0.0259189282231727  },
        { 1.14984691288476, -1.86148610164711, 0.282410704404651, -0.0137482001748907  },
        { 8.87581132716677, -6.06418672877366, 0.934174708128712, -0.047710831271599  }
    },
    // Fit for context:  AG
    {
        { -3.5155637685456, 0.445157302722387, -0.0895426643665516, 0.00489266836530986  },
        { -2.13069399948839, -0.350350414444882, 0.00136462130313272, 0.000872314444289463  },
        { 8.18160747215397, -4.75181353430953, 0.663320066629866, -0.0312212296436783  }
    },
    // Fit for context:  AT
    {
        { -0.791966407223805, -1.1076552938185, 0.121601533987201, -0.00531175337455805  },
        { -0.881444261440707, -0.731399892071115, 0.05057749339712, 9.31173789284391e-05  },
        { 1.41263899033538, -2.13869979708684, 0.29554287522677, -0.0139759898505845  }
    },
    // Fit for context:  CA
    {
        { -3.26897975135237, 0.486468075687521, -0.0760370712635174, 0.00388977870825071  },
        { -1.35534562252107, -0.638577072165448, 0.0728422429297584, -0.0017528458186216  },
        { 1.44591127557233, -1.35648438299165, 0.109211078325847, -0.00258389758690641  }
    },
    // Fit for context:  CC
    {
        { -3.06858725937395, 0.180650646816408, -0.0144155236738455, -0.000689002527617641  },
        { -6.84579191173796, 2.25516178646427, -0.339242808398693, 0.0163235455182541  },
        { 0.960878124005925, -1.79437153365345, 0.26563747809023, -0.0132138048710256  }
    },
    // Fit for context:  CG
    {
        { -2.94152963474479, 0.229999854263155, -0.0481872379951125, 0.00285101023017463  },
        { 0.869169888335733, -1.64839041486039, 0.219451593119437, -0.00974508815102731  },
        { 3.06149815060729, -2.70884449310353, 0.353919865598604, -0.0152760395255156  }
    },
    // Fit for context:  CT
    {
        { 0.713137277702561, -1.70038938696478, 0.195497764255509, -0.006277117941137  },
        { 0.547238286104262, -1.18712111898372, 0.10016377902078, 0.000101142585729921  },
        { 12.8142966669109, -7.60356150884139, 1.14279549736876, -0.0556464549926148  }
    },
    // Fit for context:  GA
    {
        { -4.04796366671756, 0.995050162374295, -0.176460865193268, 0.0098125020081521  },
        { -2.84979586974211, -0.0156903140508224, -0.0289509861158852, 0.00249431280742574  },
        { 3.58059555038287, -2.82647199241841, 0.361847788449807, -0.0165539653343998  }
    },
    // Fit for context:  GC
    {
        { -0.349844213146088, -1.43550809248254, 0.203825919576521, -0.00916693940279327  },
        { 1.63175479730304, -2.09434440700094, 0.344467568084574, -0.0177539882113115  },
        { -0.575918177053267, -1.5981774767604, 0.220590706766665, -0.0105341908671507  }
    },
    // Fit for context:  GG
    {
        { -2.7885490831891, 0.0219683824222214, -0.0434630049972039, 0.00405333456781744  },
        { 1.88946633862993, -2.07810257475806, 0.259413174779754, -0.0120217530655737  },
        { 1.69364885624261, -1.79573141305228, 0.241169169566732, -0.0114709227454174  }
    },
    // Fit for context:  GT
    {
        { -3.5677104673935, 0.358289552435854, -0.0789736924123636, 0.0036044786219107  },
        { -6.64088518072203, 2.07297640845549, -0.400256333112317, 0.0235022424593019  },
        { 4.80925899732482, -3.74080441862834, 0.510449289512398, -0.0219234011000643  }
    },
    // Fit for context:  TA
    {
        { -2.64069119760436, 0.38549684716509, -0.0816358329471849, 0.00533873493972419  },
        { 1.9113212271993, -2.61336018392704, 0.462996555125255, -0.0271165824050175  },
        { 0.238605182625403, -0.457118176718339, -0.0792840976069514, 0.00934498149921929  }
    },
    // Fit for context:  TC
    {
        { -0.392997656904618, -0.971016641524484, 0.103185513110187, -0.00226499422687784  },
        { -5.41259348922946, 1.37988834568899, -0.223433285481649, 0.0121196778796477  },
        { 0.756454968716176, -2.18832755353204, 0.294407264252205, -0.0125416902322882  }
    },
    // Fit for context:  TG
    {
        { -3.13537076489306, 0.356196913405038, -0.0801078852726234, 0.00556755944444425  },
        { 0.0417901722734827, -1.44128792816831, 0.193288926462961, -0.00976454902739811  },
        { 2.9482632980712, -2.38567585927352, 0.272478951901359, -0.0102457478965403  }
    },
    // Fit for context:  TT
    { 
        { -2.23386221240998, -0.150165813792347, 0.0503325945637299, -0.00444264897626409  },
        { 1.66643125707819, -2.11840547506457, 0.345196264109573, -0.0190747089705634  },
        { -1.9903454046927, -0.489297962270134, 0.0503913801733079, -0.00174698908180932  } 
    }};
    
inline double CalculateExpectedLogLikelihoodOfOutcomeRow(const int index, const uint8_t row, const bool secondMoment)  {
    double expectedLL = 0;
    for(size_t i = 0; i < OUTCOME_NUMBER; i++) {
        double curProb = emissionPmf[index][row][i];
        double lgCurProb = std::log(curProb);
        if(!secondMoment) {
            expectedLL +=  curProb * lgCurProb;
        } else {
            expectedLL += curProb * pow(lgCurProb, 2.0);
        }
    }
    return expectedLL;
}
    
double SP1C1PwModel::ExpectedLogLikelihoodOfOutcomeRow(const int index, const uint8_t prev, const uint8_t curr, const bool secondMoment) const {
    const auto row = (prev << 2) | curr;
    const auto moment = secondMoment ? 1 : 0;
    return cachedEmissionExpectations_[row][index][moment];
}

SP1C1PwModel::SP1C1PwModel(const SNR& snr) : snr_(snr)
{
    // Generate cached transistion probabilities
    
    double snr1 = snr_.A;
    // Adjust as necessary to clamp to bounds of SNRs observe during training
    if (snr1 < 4.0) {snr1 = 4.0;}
    else if (snr1 > 10.65) {snr1 = 10.65;}
    
    const double snr2 = snr1 * snr1, snr3 = snr2 * snr1;
    for (int ctx = 0; ctx < CONTEXT_NUMBER; ++ctx) {
        double sum = 1.0;
        ctxTrans_[ctx][0] = 1.0;
        for (size_t j = 0; j < 3; ++j) {
            double xb = transProbs[ctx][j][0] + snr1 * transProbs[ctx][j][1] +
                        snr2 * transProbs[ctx][j][2] + snr3 * transProbs[ctx][j][3];
            xb = std::exp(xb);
            ctxTrans_[ctx][j + 1] = xb;
            sum += xb;
        }
        for (size_t j = 0; j < 4; ++j)
            ctxTrans_[ctx][j] /= sum;
    }
    
    // Generate cached emission expectations
    // TODO: These are identical for all instances, either we should enrich the model or avoid doing this in a context dependent way
    for(int ctx = 0; ctx < CONTEXT_NUMBER; ctx++) {
        for (int index = 0; index < 3; index++) {
            cachedEmissionExpectations_[ctx][index][0] = CalculateExpectedLogLikelihoodOfOutcomeRow(index, ctx, false);
            cachedEmissionExpectations_[ctx][index][1] = CalculateExpectedLogLikelihoodOfOutcomeRow(index, ctx, true);
        }
    }
}

std::unique_ptr<AbstractRecursor> SP1C1PwModel::CreateRecursor(
    std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr, double scoreDiff) const
{
    return std::unique_ptr<AbstractRecursor>(
        new SP1C1PwRecursor(std::forward<std::unique_ptr<AbstractTemplate>>(tpl), mr, scoreDiff));
}

std::vector<TemplatePosition> SP1C1PwModel::Populate(const std::string& tpl) const
{
    std::vector<TemplatePosition> result;

    if (tpl.empty()) return result;

    result.reserve(tpl.size());

    // Calculate probabilities in all 16 Contexts
    uint8_t prev = detail::TranslationTable[static_cast<uint8_t>(tpl[0])];
    if (prev > 3) throw std::invalid_argument("invalid character in template!");

    for (size_t i = 1; i < tpl.size(); ++i) {
        const uint8_t curr = detail::TranslationTable[static_cast<uint8_t>(tpl[i])];
        if (curr > 3) throw std::invalid_argument("invalid character in template!");
        const auto row = (prev << 2) | curr;
        const auto params = ctxTrans_[row];
        result.emplace_back(TemplatePosition{
            tpl[i - 1], prev,
            params[0],  // match
            params[1],  // branch
            params[2],  // stick
            params[3]   // deletion
        });
        prev = curr;
    }
    result.emplace_back(TemplatePosition{tpl.back(), prev, 1.0, 0.0, 0.0, 0.0});

    return result;
}
    


double SP1C1PwModel::ExpectedLogLikelihoodForMatchEmission(uint8_t prev, uint8_t curr, bool secondMoment) const {
    return ExpectedLogLikelihoodOfOutcomeRow(static_cast<uint8_t>(MoveType::MATCH), prev, curr, secondMoment);
}
double SP1C1PwModel::ExpectedLogLikelihoodForStickEmission(uint8_t prev, uint8_t curr, bool secondMoment) const {
    return ExpectedLogLikelihoodOfOutcomeRow(static_cast<uint8_t>(MoveType::STICK), prev, curr, secondMoment);
}
double SP1C1PwModel::ExpectedLogLikelihoodForBranchEmission(uint8_t prev, uint8_t curr, bool secondMoment) const {
    return ExpectedLogLikelihoodOfOutcomeRow(static_cast<uint8_t>(MoveType::BRANCH), prev, curr, secondMoment);
}

SP1C1PwRecursor::SP1C1PwRecursor(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr,
                                 double scoreDiff)
    : Recursor<SP1C1PwRecursor>(std::forward<std::unique_ptr<AbstractTemplate>>(tpl), mr, scoreDiff)
{
}

std::vector<uint8_t> SP1C1PwRecursor::EncodeRead(const MappedRead& read)
{
    std::vector<uint8_t> result;

    result.reserve(read.Length());

    for (size_t i = 0; i < read.Length(); ++i) {
        if (read.PulseWidth[i] < 1U) throw std::runtime_error("invalid PulseWidth in read!");
        const uint8_t pw = std::min(2, read.PulseWidth[i] - 1);
        const uint8_t bp = detail::TranslationTable[static_cast<uint8_t>(read.Seq[i])];
        if (bp > 3) throw std::invalid_argument("invalid character in read!");
        const uint8_t em = (pw << 2) | bp;
        if (em > 11) throw std::runtime_error("read encoding error!");
        result.emplace_back(em);
    }

    return result;
}

double SP1C1PwRecursor::EmissionPr(MoveType move, uint8_t emission, uint8_t prev, uint8_t curr)
{
    assert(move != MoveType::DELETION);
    const auto row = (prev << 2) | curr;
    return emissionPmf[static_cast<uint8_t>(move)][row][emission] * kCounterWeight;
}

double SP1C1PwRecursor::UndoCounterWeights(const size_t nEmissions) const
{
    return -std::log(kCounterWeight) * nEmissions;
}
}  // namespace anonymous
}  // namespace Consensus
}  // namespace PacBio
