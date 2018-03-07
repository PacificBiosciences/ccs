// Author: Lance Hepler

#include "../UnanimityInternalConfig.h"

#include <cassert>
#include <cmath>
#include <memory>
#include <random>
#include <stdexcept>

#include <pacbio/consensus/ModelConfig.h>
#include <pacbio/data/Read.h>

#include "../ModelFactory.h"
#include "../Recursor.h"
#include "../Simulator.h"
#include "CounterWeight.h"
#include "HelperFunctions.h"

using namespace PacBio::Data;

namespace PacBio {
namespace Consensus {
namespace S_P1C1v2 {
namespace {

static constexpr const size_t CONTEXT_NUMBER = 16;
static constexpr const size_t OUTCOME_NUMBER = 12;

class S_P1C1v2_Model : public ModelConfig
{
public:
    static std::set<std::string> Chemistries() { return {"S/P1-C1.2", "S/P1-C1.3", "S/P2-C2"}; }
    static ModelForm Form() { return ModelForm::PWSNR; }
    S_P1C1v2_Model(const SNR& snr);
    std::unique_ptr<AbstractRecursor> CreateRecursor(const MappedRead& mr,
                                                     double scoreDiff) const override;
    std::vector<TemplatePosition> Populate(const std::string& tpl) const override;
    std::pair<Data::Read, std::vector<MoveType>> SimulateRead(
        std::default_random_engine* const rng, const std::string& tpl,
        const std::string& readname) const override;
    double ExpectedLLForEmission(MoveType move, const AlleleRep& prev, const AlleleRep& curr,
                                 MomentType moment) const override;

private:
    SNR snr_;
    double ctxTrans_[CONTEXT_NUMBER][4];
    double cachedEmissionExpectations_[CONTEXT_NUMBER][3][2];
};

// TODO(lhepler) comments regarding the CRTP
class S_P1C1v2_Recursor : public Recursor<S_P1C1v2_Recursor>
{
public:
    S_P1C1v2_Recursor(const MappedRead& mr, double scoreDiff, double counterWeight);
    static inline std::vector<uint8_t> EncodeRead(const MappedRead& read);
    inline double EmissionPr(MoveType move, uint8_t emission, const AlleleRep& prev,
                             const AlleleRep& curr) const;
    double UndoCounterWeights(size_t nEmissions) const override;

private:
    double counterWeight_;
    double nLgCounterWeight_;
};

static constexpr const double snrRanges[2][4] = {
    {4.000641, 5.654043, 4.000119, 5.751278},     // minimum
    {9.7992120, 17.958242, 9.7440690, 14.595276}  // maximum
};

static constexpr const double emissionPmf[3][CONTEXT_NUMBER][OUTCOME_NUMBER] = {
    {// matchPmf
     {0.050560033, 0.00144763859, 0.000163177865, 2.0867935e-05, 0.0726032734, 0.000842439402,
      0.000271423794, 8.43974963e-06, 0.869290866, 0.00413508802, 0.000521750699, 0.000118837204},
     {0.00462830888, 0.0139517191, 0.000101186697, 5.67650713e-05, 0.0160684007, 0.0406243136,
      0.000138793033, 5.62552642e-05, 0.0190043135, 0.904879566, 0.000174480994, 0.000295282359},
     {0.000175052112, 2.26892306e-05, 0.0196019025, 0.00400244294, 0.000509824298, 9.75492549e-06,
      0.0802293256, 0.00408560015, 0.000147604012, 0.000331175767, 0.872713612, 0.0181486146},
     {8.73952452e-05, 6.17770545e-05, 0.00903033597, 0.018502656, 0.000200885068, 9.75016622e-05,
      0.0351722645, 0.0517276852, 0.000352994508, 0.000385742081, 0.0644170155, 0.819947012},
     {0.0198831161, 0.000961303214, 0.000144597389, 4.67305367e-05, 0.0638168971, 0.000837850222,
      0.000558459661, 3.99576048e-05, 0.910672574, 0.00217152562, 0.000648278718, 0.000202152294},
     {0.0410273877, 0.0126381834, 5.80888114e-05, 4.64404749e-05, 0.0138647592, 0.0380978485,
      2.49507103e-05, 4.31344748e-05, 0.0155899299, 0.878541588, 1.74872929e-05, 3.0599793e-05},
     {8.32672171e-05, 2.20656741e-05, 0.00834519645, 0.00156910057, 0.000172765793, 2.1418975e-05,
      0.0356742714, 0.00212721578, 0.000622484386, 0.000230268885, 0.933599375, 0.0175160384},
     {3.86341223e-05, 6.93599568e-05, 0.0106638632, 0.0194476326, 2.80847249e-05, 6.80203472e-05,
      0.0350317359, 0.0561405145, 0.000248220527, 0.000735245812, 0.064117019, 0.813389881},
     {0.01055129, 0.000606336462, 0.000106046244, 6.17777656e-05, 0.0340550092, 0.000538619497,
      0.000266054642, 9.96327513e-06, 0.948832735, 0.00433935042, 0.000343332891, 0.000270289762},
     {0.00199646998, 0.0070224845, 6.10929446e-05, 5.77297092e-05, 0.0076668946, 0.0217966053,
      9.74216351e-05, 1.5916467e-05, 0.011502371, 0.949283502, 0.000259732895, 0.000225420187},
     {0.0363913971, 5.18879484e-05, 0.00586382213, 0.00119694109, 7.89913825e-05, 1.4820869e-05,
      0.0288807762, 0.00180293314, 0.000312950446, 0.000119941063, 0.909704767, 0.0155609746},
     {1.64569235e-05, 4.43086063e-05, 0.00453413579, 0.0102359104, 5.18293729e-05, 6.02931756e-05,
      0.0180114607, 0.0306943017, 0.000153656166, 0.000607235317, 0.0413053796, 0.894264433},
     {0.0214236803, 0.00145057122, 0.000196790896, 7.29502814e-05, 0.0617915495, 0.000820733636,
      0.00047074373, 5.43281455e-05, 0.908393448, 0.0044080379, 0.000772738613, 0.000118585931},
     {0.00254794588, 0.00714226772, 5.36718881e-05, 5.4246971e-05, 0.00790636565, 0.0227387232,
      3.97346769e-05, 4.53241014e-05, 0.0116566568, 0.947447228, 0.000175142411, 0.000173810276},
     {8.63087551e-05, 2.68932315e-05, 0.00961793262, 0.00169384557, 7.28762914e-05, 7.4200117e-06,
      0.037210596, 0.00195177916, 0.000458586001, 0.000265371391, 0.933488709, 0.0151035143},
     {0.0297329186, 0.000570699243, 0.00648572858, 0.00907544748, 0.00299260939, 0.000974535227,
      0.0194331272, 0.0272735285, 0.0264500987, 0.0267971751, 0.0666750407, 0.783525197}},
    {// branchPmf
     {0.271862342, 0.000108958984, 0.000108958984, 0.000108958984, 0.160747222, 0.000108958984,
      0.000108958984, 0.000108958984, 0.56586501, 0.000108958984, 0.000108958984, 0.000108958984},
     {0.000253399313, 0.0805474032, 0.000253399313, 0.000253399313, 0.000253399313, 0.11493698,
      0.000253399313, 0.000253399313, 0.000253399313, 0.800968026, 0.000253399313, 0.000253399313},
     {9.96116186e-05, 9.96116186e-05, 0.295921139, 9.96116186e-05, 9.96116186e-05, 9.96116186e-05,
      0.156000597, 9.96116186e-05, 9.96116186e-05, 9.96116186e-05, 0.546683702, 9.96116186e-05},
     {0.000266389696, 0.000266389696, 0.000266389696, 0.106949539, 0.000266389696, 0.000266389696,
      0.000266389696, 0.117892594, 0.000266389696, 0.000266389696, 0.000266389696, 0.771428411},
     {0.158494227, 3.22815102e-05, 3.22815102e-05, 3.22815102e-05, 0.135246183, 3.22815102e-05,
      3.22815102e-05, 3.22815102e-05, 0.705807649, 3.22815102e-05, 3.22815102e-05, 3.22815102e-05},
     {5.63990004e-05, 0.0396234254, 5.4605313e-05, 5.4605313e-05, 5.47353715e-05, 0.0555789327,
      5.4605313e-05, 5.4605313e-05, 5.64103869e-05, 0.904029439, 5.4605313e-05, 5.4605313e-05},
     {4.80816837e-05, 4.80816837e-05, 0.214977302, 4.80816837e-05, 4.80816837e-05, 4.80816837e-05,
      0.123027519, 4.80816837e-05, 4.80816837e-05, 4.80816837e-05, 0.661322035, 4.80816837e-05},
     {0.00019927387, 0.00019927387, 0.00019927387, 0.101763692, 0.00019927387, 0.00019927387,
      0.00019927387, 0.111141255, 0.00019927387, 0.00019927387, 0.00019927387, 0.784305219},
     {0.126620811, 3.77512454e-05, 3.77512454e-05, 3.77512454e-05, 0.0949253706, 3.77512454e-05,
      3.77512454e-05, 3.77512454e-05, 0.777925301, 3.77512454e-05, 3.77512454e-05, 3.77512454e-05},
     {8.34845269e-05, 0.0522147481, 8.34845269e-05, 8.34845269e-05, 8.34845269e-05, 0.0565042132,
      8.34845269e-05, 8.34845269e-05, 8.34845269e-05, 0.890112255, 8.34845269e-05, 8.34845269e-05},
     {0.000110815138, 0.000110427066, 0.329742326, 0.000110427066, 0.000113945866, 0.000110427066,
      0.140289483, 0.000110427066, 0.000117423649, 0.000110427066, 0.528411309, 0.000110427066},
     {0.000131337002, 0.000131337002, 0.000131337002, 0.0390710366, 0.000131337002, 0.000131337002,
      0.000131337002, 0.0617348823, 0.000131337002, 0.000131337002, 0.000131337002, 0.897355363},
     {0.121151783, 3.91167535e-05, 3.91167535e-05, 3.91167535e-05, 0.121341222, 3.91167535e-05,
      3.91167535e-05, 3.91167535e-05, 0.75695936, 3.91167535e-05, 3.91167535e-05, 3.91167535e-05},
     {6.4451422e-05, 0.030096713, 6.4451422e-05, 6.4451422e-05, 6.4451422e-05, 0.0520326862,
      6.4451422e-05, 6.4451422e-05, 6.4451422e-05, 0.916968281, 6.4451422e-05, 6.4451422e-05},
     {4.34265926e-05, 4.34265926e-05, 0.172479716, 4.34265926e-05, 4.34265926e-05, 4.34265926e-05,
      0.113044043, 4.34265926e-05, 4.34265926e-05, 4.34265926e-05, 0.713868268, 4.34265926e-05},
     {3.16836313e-05, 2.73694388e-05, 2.73694388e-05, 0.04721317, 2.77977828e-05, 2.73694388e-05,
      2.73694388e-05, 0.078916509, 3.46399946e-05, 2.73694388e-05, 2.73694388e-05, 0.873475136}},
    {// stickPmf
     {0.000103010778, 0.0380995742, 0.426353247, 0.0339918103, 0.000103010778, 0.0211011476,
      0.126299789, 0.0147406207, 0.000103010778, 0.111827288, 0.148961708, 0.0778007297},
     {0.142966063, 6.82261622e-05, 0.318828309, 0.0291750067, 0.0612651332, 6.82261622e-05,
      0.132226999, 0.0222891262, 0.108765546, 6.82261622e-05, 0.107092044, 0.0768459636},
     {0.435501353, 0.0281061159, 0.00022795498, 0.112842372, 0.0386560686, 0.0153571105,
      0.00022795498, 0.0531675523, 0.0885281789, 0.0680543734, 0.00022795498, 0.157963236},
     {0.277671243, 0.0199532671, 0.279545911, 9.68193392e-05, 0.0672566924, 0.0113474679,
      0.103986414, 9.68193392e-05, 0.0874665425, 0.0477403036, 0.104257603, 9.68193392e-05},
     {6.26766712e-05, 0.0248510469, 0.399492811, 0.0373275328, 6.26766712e-05, 0.00916381627,
      0.119657428, 0.0131950955, 6.26766712e-05, 0.212819044, 0.124778109, 0.0582137031},
     {0.125185057, 3.47389275e-05, 0.306433733, 0.0308523177, 0.0774246347, 3.52102425e-05,
      0.149641423, 0.0134852162, 0.129457219, 3.45901197e-05, 0.128579834, 0.0386692545},
     {0.314712978, 0.0245469889, 7.44182399e-05, 0.0329602878, 0.124756147, 0.0144931638,
      7.44182399e-05, 0.0161076306, 0.157254073, 0.222533644, 7.44182399e-05, 0.0920397407},
     {0.234807777, 0.0143254662, 0.196699446, 7.60379939e-05, 0.0813954759, 0.0105264534,
      0.0840430032, 7.60379939e-05, 0.112867722, 0.174110565, 0.0906157878, 7.60379939e-05},
     {0.000114280663, 0.0243444281, 0.457199123, 0.0381081135, 0.000114280663, 0.0162510324,
      0.124731361, 0.0167679501, 0.000114280663, 0.107039005, 0.107835264, 0.106809478},
     {0.118927574, 3.10562894e-05, 0.343802957, 0.0333522609, 0.0496397906, 3.10562894e-05,
      0.195022097, 0.0163660717, 0.0539508263, 3.10562894e-05, 0.125648942, 0.0630410307},
     {0.335673887, 0.0267972214, 0.00093244908, 0.0455030532, 0.129837105, 0.0126090082,
      0.000180646254, 0.0192047603, 0.142185037, 0.0744691665, 0.000164140614, 0.211647424},
     {0.299454906, 0.0203545133, 0.270281208, 0.000115353453, 0.10064444, 0.00705182304,
      0.114896506, 0.000115353453, 0.100554091, 0.0384270526, 0.0474126337, 0.000115353453},
     {7.92484966e-05, 0.0190180551, 0.262244947, 0.019509475, 7.92484966e-05, 0.0113606249,
      0.0949860067, 0.0183237604, 7.92484966e-05, 0.106081588, 0.21160129, 0.256240265},
     {0.125916528, 4.15946843e-05, 0.209865333, 0.0203882413, 0.0577305702, 4.15946843e-05,
      0.0904461035, 0.018649162, 0.0587107584, 4.15946843e-05, 0.129821536, 0.28813901},
     {0.33820122, 0.0210315415, 0.000104202635, 0.0409843175, 0.118326009, 0.0105217011,
      0.000104202635, 0.0250824494, 0.10336307, 0.0371536554, 0.000104202635, 0.304502416},
     {0.148349765, 0.0127413602, 0.22312407, 4.96955324e-05, 0.0595806675, 0.00463335009,
      0.177668167, 4.83714586e-05, 0.0574131588, 0.0412555016, 0.274858141, 4.85075474e-05}}};

static constexpr const double transProbs[CONTEXT_NUMBER][3][4] = {
    {// AA
     {5.85956362, -4.41463293, 0.683309462, -0.0346718671},
     {7.06499623, -4.40259797, 0.607137416, -0.027878395},
     {1.73648922, -0.935476496, 0.0264529997, 0.00196800649}},
    {// AC
     {2.68872981, -1.33398847, 0.0860003242, -0.00189734598},
     {3.55997286, -1.60728992, 0.130435583, -0.00340039547},
     {4.95711278, -1.98021021, 0.144128237, -0.00360784403}},
    {// AG
     {-3.23979278, 0.356138788, -0.0870220245, 0.00546284744},
     {2.5291947, -2.13705602, 0.206929395, -0.00511956903},
     {2.20068033, -0.996606761, -0.00278577326, 0.00492497905}},
    {// AT
     {-1.15847683, -0.771986547, 0.0645814917, -0.0018601879},
     {-2.78585757, 0.178618309, -0.0499252647, 0.00257301448},
     {3.08433285, -1.67474855, 0.14589205, -0.00469624396}},
    {// CA
     {-0.451942513, -0.575360527, 0.0401324889, 0.000694689568},
     {0.58830656, -1.53028855, 0.204160977, -0.00810273241},
     {2.35067779, -1.31686514, 0.0730142817, -1.00923273e-05}},
    {// CC
     {-7.2469072, 1.26130953, -0.111198017, 0.0031739936},
     {-5.43347752, 0.980596171, -0.0930509844, 0.00284106343},
     {-4.62125702, 0.536176189, -0.0532120198, 0.00156457856}},
    {// CG
     {-1.85080654, -0.271149887, 0.0233007441, -0.000158657488},
     {3.74845392, -2.86890561, 0.378280475, -0.0155172951},
     {5.12582688, -3.00674655, 0.342928086, -0.0137785283}},
    {// CT
     {1.59462161, -1.32139342, 0.112979516, -0.00350654749},
     {-2.66688067, 0.149169695, -0.0394165557, 0.00221529117},
     {3.35981388, -1.6812294, 0.144580211, -0.00463679475}},
    {// GA
     {1.94929663, -1.87204746, 0.262415331, -0.0114952071},
     {0.776645323, -1.57584288, 0.170612002, -0.00429563295},
     {9.14677248, -5.01761335, 0.672305833, -0.031991748}},
    {// GC
     {10.8094257, -3.57210251, 0.294531042, -0.0079448247},
     {-1.6188921, -0.212323741, 0.0145898135, -0.000189460808},
     {0.516468575, -0.980880476, 0.0607592136, -0.00125531492}},
    {// GG
     {-1.967728, -0.522486336, 0.0601164523, -0.00183881474},
     {3.12610849, -2.59191292, 0.325979211, -0.0138882827},
     {3.05491681, -2.21846076, 0.287404314, -0.0128841407}},
    {// GT
     {-5.79447669, 0.996639643, -0.120105223, 0.00438200289},
     {7.06158739, -2.93404754, 0.268804061, -0.00798584198},
     {6.70262773, -2.87567511, 0.264224627, -0.00821191872}},
    {// TA
     {-7.07222018, 2.94211219, -0.543110434, 0.0320758432},
     {-0.220904721, -0.873637208, 0.06909082, 0.000725498815},
     {3.06972721, -1.54851697, 0.0802862366, 0.00176909841}},
    {// TC
     {2.13319889, -1.20882369, 0.0954386191, -0.00244438603},
     {-2.18846256, 0.022602406, -0.0121265739, 0.000696342338},
     {-5.49462196, 0.733656731, -0.098767722, 0.00362214647}},
    {// TG
     {0.459855955, -1.41666334, 0.205580055, -0.00922652629},
     {1.40013661, -1.69751326, 0.165896712, -0.00352395246},
     {3.64491162, -2.25206115, 0.23605282, -0.00926403736}},
    {// TT
     {-4.66331239, 0.6977867, -0.0659908318, 0.00212194439},
     {-0.800912389, -0.311917271, 0.00243369751, 0.000960389041},
     {6.01995449, -2.74903881, 0.271552639, -0.00905647598}}};

inline double CalculateExpectedLLForEmission(const size_t move, const uint8_t row,
                                             const size_t moment)
{
    double expectedLL = 0;
    for (size_t i = 0; i < OUTCOME_NUMBER; i++) {
        double curProb = emissionPmf[move][row][i];
        double lgCurProb = std::log(curProb);
        if (moment == static_cast<uint8_t>(MomentType::FIRST))
            expectedLL += curProb * lgCurProb;
        else if (moment == static_cast<uint8_t>(MomentType::SECOND))
            expectedLL += curProb * (lgCurProb * lgCurProb);
    }
    return expectedLL;
}

S_P1C1v2_Model::S_P1C1v2_Model(const SNR& snr)
    : snr_(ClampSNR(snr, SNR{snrRanges[0]}, SNR{snrRanges[1]}))
{
    for (size_t ctx = 0; ctx < CONTEXT_NUMBER; ++ctx) {
        const uint8_t bp = ctx & 3;  // (equivalent to % 4)
        const double snr1 = snr_[bp], snr2 = snr1 * snr1, snr3 = snr2 * snr1;
        double sum = 1.0;

        // cached transitions
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

        // cached expectations
        for (size_t move = 0; move < 3; ++move)
            for (size_t moment = 0; moment < 2; ++moment)
                cachedEmissionExpectations_[ctx][move][moment] =
                    CalculateExpectedLLForEmission(move, ctx, moment);
    }
}

std::unique_ptr<AbstractRecursor> S_P1C1v2_Model::CreateRecursor(const MappedRead& mr,
                                                                 double scoreDiff) const
{
    const double counterWeight = CounterWeight(
        [this](size_t ctx, MoveType m) { return ctxTrans_[ctx][static_cast<uint8_t>(m)]; },
        [](size_t ctx, MoveType m) {
            double r = 0.0;
            for (size_t o = 0; o < OUTCOME_NUMBER; ++o) {
                const double p = emissionPmf[static_cast<uint8_t>(m)][ctx][o];
                if (p > 0.0) r += p * std::log(p);
            }
            return r;
        },
        CONTEXT_NUMBER);

    return std::make_unique<S_P1C1v2_Recursor>(mr, scoreDiff, counterWeight);
}

std::vector<TemplatePosition> S_P1C1v2_Model::Populate(const std::string& tpl) const
{
    auto rowFetcher = [this](const NCBI2na prev, const NCBI2na curr) -> const double(&)[4]
    {
        const auto row = EncodeContext16(prev, curr);
        const double(&params)[4] = ctxTrans_[row];
        return params;
    };
    return AbstractPopulater(tpl, rowFetcher);
}

double S_P1C1v2_Model::ExpectedLLForEmission(const MoveType move, const AlleleRep& prev,
                                             const AlleleRep& curr, const MomentType moment) const
{
    auto cachedEmissionVisitor = [this](const MoveType move, const NCBI2na prev, const NCBI2na curr,
                                        const MomentType moment) -> double {
        const auto row = EncodeContext16(prev, curr);
        return cachedEmissionExpectations_[row][static_cast<uint8_t>(move)]
                                          [static_cast<uint8_t>(moment)];
    };
    return AbstractExpectedLLForEmission(move, prev, curr, moment, cachedEmissionVisitor);
}

S_P1C1v2_Recursor::S_P1C1v2_Recursor(const MappedRead& mr, double scoreDiff, double counterWeight)
    : Recursor<S_P1C1v2_Recursor>(mr, scoreDiff)
    , counterWeight_{counterWeight}
    , nLgCounterWeight_{-std::log(counterWeight_)}
{
}

std::vector<uint8_t> S_P1C1v2_Recursor::EncodeRead(const MappedRead& read)
{
    std::vector<uint8_t> result;
    result.reserve(read.Length());

    for (size_t i = 0; i < read.Length(); ++i) {
        result.emplace_back(EncodeBase(read.Seq[i], read.PulseWidth[i]));
    }

    return result;
}

double S_P1C1v2_Recursor::EmissionPr(const MoveType move, const uint8_t emission,
                                     const AlleleRep& prev, const AlleleRep& curr) const
{
    return AbstractEmissionPr(emissionPmf, move, emission, prev, curr) * counterWeight_;
}

double S_P1C1v2_Recursor::UndoCounterWeights(const size_t nEmissions) const
{
    return nLgCounterWeight_ * nEmissions;
}

inline std::pair<Data::SNR, std::vector<TemplatePosition>> S_P1C1v2_InitialiseModel(
    std::default_random_engine* const rng, const std::string& tpl)
{
    Data::SNR snrs{0, 0, 0, 0};
    for (uint8_t i = 0; i < 4; ++i) {
        snrs[i] = std::uniform_real_distribution<double>{snrRanges[0][i], snrRanges[1][i]}(*rng);
    }

    const S_P1C1v2_Model model{snrs};
    std::vector<TemplatePosition> transModel = model.Populate(tpl);

    return {snrs, transModel};
}

BaseData S_P1C1v2_GenerateReadData(std::default_random_engine* const rng, const MoveType state,
                                   const AlleleRep& prev, const AlleleRep& curr)
{
    // distribution is arbitrary at the moment, as
    // IPD is not a covariate of the consensus HMM
    std::uniform_int_distribution<uint8_t> ipdDistrib(1, 5);

    std::array<double, OUTCOME_NUMBER> emissionDist;
    for (size_t i = 0; i < OUTCOME_NUMBER; ++i) {
        emissionDist[i] = AbstractEmissionPr(emissionPmf, state, i, prev, curr);
    }

    std::discrete_distribution<uint8_t> outcomeDistrib(emissionDist.cbegin(), emissionDist.cend());

    const uint8_t event = outcomeDistrib(*rng);
    const std::pair<char, uint8_t> outcome = DecodeEmission(event);

    return {outcome.first, outcome.second, ipdDistrib(*rng)};
}

std::pair<Data::Read, std::vector<MoveType>> S_P1C1v2_Model::SimulateRead(
    std::default_random_engine* const rng, const std::string& tpl,
    const std::string& readname) const
{
    return SimulateReadImpl(rng, tpl, readname, S_P1C1v2_InitialiseModel,
                            S_P1C1v2_GenerateReadData);
}

}  // namespace anonymous
}  // namespace S_P1C1v2

REGISTER_MODEL_IMPL(S_P1C1v2)

}  // namespace Consensus
}  // namespace PacBio
