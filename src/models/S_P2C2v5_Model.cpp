// Author: Lance Hepler

#include "../UnanimityInternalConfig.h"

#include <cassert>
#include <cmath>
#include <memory>
#include <random>
#include <stdexcept>

#include <pacbio/consensus/ModelConfig.h>
#include <pacbio/data/Read.h>
#include <pacbio/data/internal/BaseEncoding.h>

#include "../ModelFactory.h"
#include "../Recursor.h"
#include "../Simulator.h"
#include "CounterWeight.h"
#include "HelperFunctions.h"

using namespace PacBio::Data;

namespace PacBio {
namespace Consensus {
namespace S_P2C2v5 {
namespace {

static constexpr const size_t CONTEXT_NUMBER = 16;
static constexpr const size_t OUTCOME_NUMBER = 12;

class S_P2C2v5_Model : public ModelConfig
{
public:
    static std::set<std::string> Chemistries() { return {"S/P2-C2/5.0"}; }
    static ModelForm Form() { return ModelForm::PWSNR; }
    S_P2C2v5_Model(const SNR& snr);
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
class S_P2C2v5_Recursor : public Recursor<S_P2C2v5_Recursor>
{
public:
    S_P2C2v5_Recursor(const MappedRead& mr, double scoreDiff, double counterWeight);
    static inline std::vector<uint8_t> EncodeRead(const MappedRead& read);
    inline double EmissionPr(MoveType move, uint8_t emission, const AlleleRep& prev,
                             const AlleleRep& curr) const;
    double UndoCounterWeights(size_t nEmissions) const override;

private:
    double counterWeight_;
    double nLgCounterWeight_;
};

static constexpr const double snrRanges[2][4] = {
    {3.91070819, 7.37540293, 3.84641552, 6.27542830},  // minimum
    {8.06702614, 15.2471046, 7.02776623, 11.3469543}   // maximum
};

static constexpr const double emissionPmf[3][CONTEXT_NUMBER][OUTCOME_NUMBER] = {
    {// matchPmf
     {0, 0.000180313185, 0, 1.01091615e-06, 0.00352197085, 0.019010118, 7.33665722e-06,
      2.45437478e-05, 0.964232358, 0.012245689, 0.000606040406, 0.000170618772},
     {0, 0.00163157289, 0, 1.61106086e-05, 0.0002339041, 0.015752956, 1.98545947e-07,
      6.14642833e-05, 0.00292805283, 0.979067531, 0.000265869341, 4.23399988e-05},
     {0, 5.11199878e-08, 0, 8.81535441e-05, 6.09701703e-06, 0.000120041542, 0.00418106818,
      0.0186350696, 0.000143659898, 0.000337260556, 0.928654474, 0.0478341246},
     {0, 4.20774984e-06, 0, 0.00100463255, 6.20932357e-06, 0.000200913588, 0.000814896042,
      0.0211494868, 7.04763101e-05, 0.000634728089, 0.0140405685, 0.962073881},
     {0, 0.000100786916, 0, 1.06397215e-07, 0.00422094102, 0.0198224993, 1.32150878e-05,
      0.000126391143, 0.963236551, 0.0113698784, 0.000862125915, 0.000247505016},
     {0, 0.00150409409, 0, 2.24863787e-05, 0.000219495725, 0.0168965224, 1.5963561e-07,
      9.7206503e-07, 0.00170332286, 0.979650331, 1.64853182e-06, 9.67668332e-07},
     {0, 5.40174236e-06, 0, 1.2289404e-05, 2.51090942e-06, 0.00013358247, 0.00183814326,
      0.00846971697, 0.000867641038, 0.000299484643, 0.952805673, 0.0355655567},
     {0, 4.0944742e-06, 0, 0.00144215522, 2.1456507e-05, 0.000143220285, 0.000801561351,
      0.0256747187, 0.000320224499, 0.000404895805, 0.0152899374, 0.955897736},
     {0, 7.45953044e-05, 0, 1.88483032e-11, 0.00183129742, 0.00952630737, 3.79438695e-06,
      8.41250377e-05, 0.978487828, 0.00963010363, 5.17630549e-05, 0.000310185962},
     {0, 0.000707729559, 0, 3.57413036e-07, 0.000184684294, 0.00877226678, 8.4738145e-06,
      4.8609979e-06, 0.00204433176, 0.987745911, 0.000377108775, 0.000154275681},
     {0, 8.94639403e-06, 0, 6.76553774e-05, 9.08896585e-06, 9.2579371e-06, 0.0018172668,
      0.00743661276, 0.000590718639, 0.000200250166, 0.958113701, 0.0317465019},
     {0, 1.10286805e-08, 0, 0.00051420915, 2.73045935e-08, 0.000150944211, 0.000475102224,
      0.0120964654, 0.000156045677, 0.000743346818, 0.00871952412, 0.977144324},
     {0, 0.000106283202, 0, 6.14989186e-06, 0.00297849039, 0.0157756968, 5.93862254e-07,
      6.24672993e-05, 0.967332481, 0.0129764955, 0.00071307378, 4.82679864e-05},
     {0, 0.000711743788, 0, 4.29817416e-12, 0.000138858898, 0.00834332455, 1.07763728e-07,
      3.98076484e-07, 0.00183089129, 0.988776882, 0.000148575774, 4.92174644e-05},
     {0, 1.47749994e-06, 0, 3.74888487e-05, 1.50070824e-05, 1.31657347e-05, 0.0020042872,
      0.00824849623, 0.00104092736, 0.000381194952, 0.956556364, 0.0317015912},
     {0, 1.32369096e-06, 0, 0.000579281915, 6.77425403e-07, 1.03383747e-05, 0.000353223056,
      0.0107239563, 1.23813368e-06, 8.81159604e-06, 0.00831039342, 0.980010756}},
    {// branchPmf
     {0, 0, 0, 0, 0.0761051598, 0, 0, 0, 0.92389484, 0, 0, 0},
     {0, 0.0133502664, 0, 0, 0, 0.137719726, 0, 0, 0, 0.848930007, 0, 0},
     {0, 0, 0, 0, 0, 0, 0.0203164781, 0, 0, 0, 0.979683522, 0},
     {0, 0, 0, 0.00956062646, 0, 0, 0, 0.130634211, 0, 0, 0, 0.859805162},
     {0, 0, 0, 0, 0.0120662791, 0, 0, 0, 0.987933721, 0, 0, 0},
     {0, 0.00754649446, 0, 0, 0, 0.0535200842, 0, 0, 0, 0.938933421, 0, 0},
     {0, 0, 0, 0, 0, 0, 0.0139128638, 0, 0, 0, 0.986087136, 0},
     {0, 0, 0, 0.00815306979, 0, 0, 0, 0.0971599384, 0, 0, 0, 0.894686992},
     {0, 0, 0, 0, 0.00842044346, 0, 0, 0, 0.991579557, 0, 0, 0},
     {0, 0.0101739992, 0, 0, 0, 0.0742441474, 0, 0, 0, 0.915581853, 0, 0},
     {0, 0, 0, 0, 0, 0, 0.0438224581, 0, 0, 0, 0.956177542, 0},
     {0, 0, 0, 0.00548230676, 0, 0, 0, 0.0610015875, 0, 0, 0, 0.933516106},
     {0, 0, 0, 0, 0.007524987, 0, 0, 0, 0.992475013, 0, 0, 0},
     {0, 0.00526638542, 0, 0, 0, 0.047146029, 0, 0, 0, 0.947587586, 0, 0},
     {0, 0, 0, 0, 0, 0, 0.00990403989, 0, 0, 0, 0.99009596, 0},
     {0, 0, 0, 0.00285972488, 0, 0, 0, 0.044549056, 0, 0, 0, 0.952591219}},
    {// stickPmf
     {0, 0.00627065742, 0, 0.00577886791, 0, 0.0961851139, 0.0196517555, 0.0623361085, 0,
      0.300764084, 0.342927482, 0.16608593},
     {0, 0, 0, 0.00466123301, 0.00758966825, 0, 0.0369524672, 0.117726338, 0.171519161, 0,
      0.419326418, 0.242224714},
     {0, 0.0162737951, 0, 0.0180673462, 0.00344712331, 0.0848166505, 0, 0.178427923, 0.00176079705,
      0.164047437, 0, 0.533158929},
     {0, 0.0183109251, 0, 0, 0.0191880875, 0.145505518, 0.0353251477, 0, 0.00215542401, 0.23830516,
      0.541209738, 0},
     {0, 0.00471013932, 0, 0.00622377965, 0, 0.136610846, 0.0162629592, 0.051051271, 0, 0.42112594,
      0.241957891, 0.122057174},
     {0, 0, 0, 0.0040276349, 0.00882937099, 0, 0.0271793145, 0.0898693705, 0.279788063, 0,
      0.409646841, 0.180659405},
     {0, 0.00564989198, 0, 0.00777629437, 0.0200504759, 0.06581245, 0, 0.126974941, 0.346675962,
      0.0799152957, 0, 0.34714469},
     {0, 0.00901392552, 0, 0, 0.0156566886, 0.0692888812, 0.0160658049, 0, 0.285234034, 0.329661799,
      0.275078867, 0},
     {0, 0.00762020132, 0, 0.00916426472, 0, 0.202689918, 0.0208904777, 0.104985428, 0, 0.277883402,
      0.104218462, 0.272547847},
     {0, 0, 0, 0.0060586789, 0.00840759497, 0, 0.0392703088, 0.164083818, 0.151353568, 0,
      0.350816456, 0.280009576},
     {0, 0.00594188465, 0, 0.00387367678, 0.0170527062, 0.0539388568, 0, 0.0678775793, 0.267335382,
      0.134832445, 0, 0.449147469},
     {0, 0.0210910872, 0, 0, 0.0343929347, 0.14510557, 0.0359179443, 0, 0.533553619, 0.213513461,
      0.0164253844, 0},
     {0, 0.00545163945, 0, 0.00299665655, 0, 0.114241543, 0.00920617524, 0.0235287971, 0,
      0.257362755, 0.234305451, 0.352906982},
     {0, 0, 0, 0.00333771279, 0.00894873436, 0, 0.0186910645, 0.052004696, 0.170146813, 0,
      0.272438848, 0.474432131},
     {0, 0.00442107785, 0, 0.00724239931, 0.017066735, 0.0663295994, 0, 0.10731725, 0.324097388,
      0.103699271, 0, 0.369826279},
     {0, 0.0059224117, 0, 0, 0.00934985504, 0.0356605899, 0.0149281548, 0, 0.269045679, 0.31090005,
      0.35419326, 0}}};

static constexpr const double transProbs[CONTEXT_NUMBER][3][4] = {
    {// AA
     {-4.59501616, -0.703891003, 0.13968334, -0.00910714009},
     {-0.835790459, -1.17248459, 0.124455964, -0.00401216571},
     {3.88047539, -2.33656165, 0.290541278, -0.011987587}},
    {// AC
     {-2.5356158, -0.190661406, -0.00315560034, 0.000680398602},
     {-5.11602184, 0.444105013, -0.0467995812, 0.00158711075},
     {-3.20478858, -0.0507252816, -0.00694193567, 0.000524056995}},
    {// AG
     {-2.95775917, 0.0434562694, -0.051753628, 0.00465543512},
     {1.34265907, -1.93402064, 0.175294352, -0.0027619221},
     {4.19947351, -2.64955863, 0.335166047, -0.0141807296}},
    {// AT
     {-3.88326282, 0.147445866, -0.0470401835, 0.00255512945},
     {4.22276046, -2.57671776, 0.268217773, -0.0101211323},
     {6.69686611, -2.95713855, 0.282464375, -0.00887025284}},
    {// CA
     {-2.5337278, -0.00694836499, -0.00304256722, 0.000723936549},
     {-3.06424651, -0.212172869, 0.0369297893, -0.00167571305},
     {4.76782904, -2.90827139, 0.365452837, -0.0149534154}},
    {// CC
     {-5.09916013, 0.469087647, -0.0281528727, 0.000554733756},
     {-5.33795356, 0.562756071, -0.0371942269, 0.000737215698},
     {-0.243107121, -0.522616567, 0.0365308646, -0.000838424089}},
    {// CG
     {-1.34013506, -0.792612921, 0.12372833, -0.00710691444},
     {2.69256994, -2.87110628, 0.437346106, -0.0227573276},
     {2.77691143, -2.31227632, 0.282607776, -0.0108952596}},
    {// CT
     {-2.68716866, 0.043111271, -0.0340305562, 0.00201357788},
     {-0.748181286, -0.657754046, 0.0489520693, -0.00129953423},
     {2.57646944, -1.48666035, 0.119489274, -0.00292702644}},
    {// GA
     {-1.65568797, -0.570179253, 0.100667528, -0.00550569443},
     {-1.61544634, -0.91563415, 0.109431411, -0.00416552791},
     {1.68408285, -1.72829129, 0.178884968, -0.00512357326}},
    {// GC
     {-2.3424541, -0.167735889, 0.00581487311, 6.21689224e-05},
     {-2.33896021, -0.205799672, 0.0111451399, -7.63565098e-05},
     {-1.21421339, -0.656491402, 0.0448956402, -0.000872934449}},
    {// GG
     {-7.10724415, 0.550738776, 0.00210501067, -0.00370835465},
     {-2.55568624, 0.232223355, -0.163503979, 0.0134946694},
     {1.12876345, -1.18404958, 0.112125326, -0.00202495079}},
    {// GT
     {2.91074017, -1.84710464, 0.175699499, -0.0054671922},
     {0.862353407, -1.25242227, 0.07860008, -0.00095188595},
     {1.60740721, -1.23052584, 0.0826531011, -0.00125966513}},
    {// TA
     {0.82120627, -1.56434011, 0.249120429, -0.0122513964},
     {0.746784419, -1.8647586, 0.27924597, -0.0122946814},
     {6.35622429, -3.83920958, 0.514244254, -0.0225187565}},
    {// TC
     {-3.69160997, 0.198856011, -0.0176361137, 0.000556398196},
     {-4.86462797, 0.33853425, -0.0284353557, 0.000863706132},
     {0.615283617, -1.16555118, 0.0920056813, -0.00230287836}},
    {// TG
     {-0.895012204, -0.671731478, 0.0504996711, 0.000942470722},
     {-0.545390091, -1.2670154, 0.129085712, -0.00329057752},
     {2.11001811, -1.66049057, 0.127729869, -0.000471511975}},
    {// TT
     {2.80777791, -1.69864481, 0.202000941, -0.00776406841},
     {3.15779542, -1.82508282, 0.176390223, -0.00572475801},
     {-0.0992454986, -0.601417099, 0.033289302, -0.000168035002}}};

inline double CalculateExpectedLLForEmission(const size_t move, const uint8_t row,
                                             const size_t moment)
{
    double expectedLL = 0;
    for (size_t i = 0; i < OUTCOME_NUMBER; i++) {
        double curProb = emissionPmf[move][row][i];
        double lgCurProb = std::log(curProb);
        if (!std::isfinite(lgCurProb)) continue;
        if (moment == static_cast<uint8_t>(MomentType::FIRST))
            expectedLL += curProb * lgCurProb;
        else if (moment == static_cast<uint8_t>(MomentType::SECOND))
            expectedLL += curProb * (lgCurProb * lgCurProb);
    }
    return expectedLL;
}

S_P2C2v5_Model::S_P2C2v5_Model(const SNR& snr)
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

std::unique_ptr<AbstractRecursor> S_P2C2v5_Model::CreateRecursor(const MappedRead& mr,
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

    return std::make_unique<S_P2C2v5_Recursor>(mr, scoreDiff, counterWeight);
}

std::vector<TemplatePosition> S_P2C2v5_Model::Populate(const std::string& tpl) const
{
    auto rowFetcher = [this](const NCBI2na prev, const NCBI2na curr) -> const double(&)[4]
    {
        const auto row = EncodeContext16(prev, curr);
        const double(&params)[4] = ctxTrans_[row];
        return params;
    };
    return AbstractPopulater(tpl, rowFetcher);
}

double S_P2C2v5_Model::ExpectedLLForEmission(const MoveType move, const AlleleRep& prev,
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

S_P2C2v5_Recursor::S_P2C2v5_Recursor(const MappedRead& mr, double scoreDiff, double counterWeight)
    : Recursor<S_P2C2v5_Recursor>(mr, scoreDiff)
    , counterWeight_{counterWeight}
    , nLgCounterWeight_{-std::log(counterWeight_)}
{
}

std::vector<uint8_t> S_P2C2v5_Recursor::EncodeRead(const MappedRead& read)
{
    std::vector<uint8_t> result;
    result.reserve(read.Length());

    for (size_t i = 0; i < read.Length(); ++i) {
        result.emplace_back(EncodeBase(read.Seq[i], read.PulseWidth[i]));
    }

    return result;
}

double S_P2C2v5_Recursor::EmissionPr(const MoveType move, const uint8_t emission,
                                     const AlleleRep& prev, const AlleleRep& curr) const
{
    return AbstractEmissionPr(emissionPmf, move, emission, prev, curr) * counterWeight_;
}

double S_P2C2v5_Recursor::UndoCounterWeights(const size_t nEmissions) const
{
    return nLgCounterWeight_ * nEmissions;
}

inline std::pair<Data::SNR, std::vector<TemplatePosition>> S_P2C2v5_InitialiseModel(
    std::default_random_engine* const rng, const std::string& tpl)
{
    Data::SNR snrs{0, 0, 0, 0};
    for (uint8_t i = 0; i < 4; ++i) {
        snrs[i] = std::uniform_real_distribution<double>{snrRanges[0][i], snrRanges[1][i]}(*rng);
    }

    const S_P2C2v5_Model model{snrs};
    std::vector<TemplatePosition> transModel = model.Populate(tpl);

    return {snrs, transModel};
}

BaseData S_P2C2v5_GenerateReadData(std::default_random_engine* const rng, const MoveType state,
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

std::pair<Data::Read, std::vector<MoveType>> S_P2C2v5_Model::SimulateRead(
    std::default_random_engine* const rng, const std::string& tpl,
    const std::string& readname) const
{
    return SimulateReadImpl(rng, tpl, readname, S_P2C2v5_InitialiseModel,
                            S_P2C2v5_GenerateReadData);
}

}  // namespace anonymous
}  // namespace S_P2C2v5

REGISTER_MODEL_IMPL(S_P2C2v5)

}  // namespace Consensus
}  // namespace PacBio
