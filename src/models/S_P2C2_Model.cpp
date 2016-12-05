// Copyright (c) 2014-2016, Pacific Biosciences of California, Inc.
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

// Author: Lance Hepler

#include <cassert>
#include <cmath>
#include <memory>
#include <stdexcept>

#include <pacbio/consensus/ModelConfig.h>
#include <pacbio/data/Read.h>

#include "../ModelFactory.h"
#include "../Recursor.h"

using namespace PacBio::Data;

namespace PacBio {
namespace Consensus {
namespace {

template <typename T>
inline T clip(const T val, const T range[2])
{
    return std::max(range[0], std::min(val, range[1]));
}

constexpr double kCounterWeight = 20.0;
constexpr size_t nOutcomes = 12;
constexpr size_t nContexts = 16;

class S_P2C2_Model : public ModelConfig
{
    REGISTER_MODEL(S_P2C2_Model);

public:
    static std::set<std::string> Names() { return {"S/P2-C2::PwSnr"}; }
    S_P2C2_Model(const SNR& snr);
    std::unique_ptr<AbstractRecursor> CreateRecursor(std::unique_ptr<AbstractTemplate>&& tpl,
                                                     const MappedRead& mr, double scoreDiff) const;
    std::vector<TemplatePosition> Populate(const std::string& tpl) const;
    double ExpectedLLForEmission(MoveType move, uint8_t prev, uint8_t curr,
                                 MomentType moment) const;

private:
    SNR snr_;
    double ctxTrans_[nContexts][4];
    double cachedEmissionExpectations_[nContexts][3][2];
};

REGISTER_MODEL_IMPL(S_P2C2_Model);

// TODO(lhepler) comments regarding the CRTP
class S_P2C2_Recursor : public Recursor<S_P2C2_Recursor>
{
public:
    S_P2C2_Recursor(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr,
                    double scoreDiff);
    static inline std::vector<uint8_t> EncodeRead(const MappedRead& read);
    static inline double EmissionPr(MoveType move, uint8_t emission, uint8_t prev, uint8_t curr);
    virtual double UndoCounterWeights(size_t nEmissions) const;
};

constexpr double snrRanges[4][2] = {{4.06761408, 7.87333488},   // A
                                    {7.56100464, 14.9044018},   // C
                                    {4.05258036, 6.67833996},   // G
                                    {6.47991323, 10.7843523}};  // T

constexpr double emissionPmf[3][nContexts][nOutcomes] = {
    {// matchPmf
     {0.00735664387, 0.00190428659, 0.000128987185, 3.19609285e-05, 0.0579789127, 0.00159365177,
      0.000377851288, 2.64641792e-05, 0.925861388, 0.0039544046, 0.000740333577, 4.51157354e-05},
     {0.00268702382, 0.0103554239, 5.5569919e-05, 0.00010710572, 0.0119429479, 0.0349107323,
      0.000173873688, 3.71404194e-05, 0.0132184981, 0.925935393, 0.000163937702, 0.000412353133},
     {0.000303779953, 5.76150405e-05, 0.00682147837, 0.00396954594, 0.000590520182, 1.49898804e-05,
      0.0623840757, 0.00479826006, 0.000162992201, 0.000111189159, 0.896219743, 0.0245658107},
     {0.000142620524, 3.85177333e-05, 0.00440064256, 0.018046601, 0.000210517258, 2.23325708e-05,
      0.0272865144, 0.0538911959, 0.000157354862, 0.000143746285, 0.0491787167, 0.84648124},
     {0.00733096773, 0.00166007937, 0.000121497496, 0.000116125734, 0.0573233988, 0.00126886069,
      0.00071151983, 7.74436926e-05, 0.927257837, 0.0023603384, 0.00153974443, 0.000232186399},
     {0.00219335118, 0.0112686045, 4.01328525e-05, 5.02885136e-05, 0.0117893856, 0.0371141169,
      1.68691479e-05, 1.56068136e-05, 0.0116461848, 0.925851525, 1.25205575e-05, 1.41414695e-06},
     {0.000117348076, 1.75083299e-05, 0.00297146156, 0.0016269759, 0.000296647296, 4.87166628e-05,
      0.0284358458, 0.002344226, 0.000568544223, 0.000237602686, 0.942157699, 0.021177424},
     {9.37972006e-05, 6.24310526e-05, 0.00526664739, 0.0207263945, 0.00012215642, 2.11772913e-05,
      0.0320582675, 0.06180119, 0.000243487269, 0.000476247908, 0.0530240115, 0.826104192},
     {0.00314264389, 0.00082820948, 9.08246005e-05, 9.38915756e-05, 0.028554939, 0.000783143672,
      0.000364350553, 2.37054476e-05, 0.961214998, 0.00421586002, 0.000391376658, 0.000296057041},
     {0.00119090385, 0.00593803619, 2.03588382e-05, 4.41352391e-05, 0.00628544188, 0.0204138621,
      0.000138312951, 6.76020532e-05, 0.00866986949, 0.956464996, 0.000533431456, 0.000233049758},
     {6.43734294e-05, 1.22650033e-05, 0.00210011572, 0.00131543919, 9.48512735e-05, 9.02245873e-06,
      0.0243938729, 0.00215317529, 0.00018372826, 3.35949921e-05, 0.950122398, 0.0195171634},
     {3.66710321e-05, 4.51849068e-05, 0.0024483006, 0.010397222, 8.29270369e-05, 1.71244431e-05,
      0.0166213186, 0.0346020058, 8.05756611e-05, 0.000443816688, 0.0317772895, 0.903447564},
     {0.00584782932, 0.00152030367, 8.92185665e-05, 8.75578495e-05, 0.0450024718, 0.00128621735,
      0.000488191126, 3.73830825e-05, 0.940618603, 0.00388969117, 0.00110426833, 2.82644734e-05},
     {0.00124169992, 0.00586665604, 3.58345575e-05, 7.76586763e-05, 0.00564348936, 0.0178975167,
      0.000131182094, 2.75294579e-05, 0.00733345492, 0.961559728, 0.000111068097, 7.41824771e-05},
     {5.36955022e-05, 3.60866503e-05, 0.00270197625, 0.00158018444, 8.86444887e-05, 7.11501847e-06,
      0.0271717889, 0.00206535023, 0.000463746178, 0.000180875439, 0.946297025, 0.0193535123},
     {4.54423247e-05, 1.30264171e-05, 0.00201775615, 0.00868651464, 2.93835524e-05, 8.6504974e-06,
      0.0146074997, 0.0295302164, 1.74359077e-05, 6.74859902e-06, 0.0306931623, 0.914344164}},
    {// branchPmf
     {0.13965274, 0.000214059166, 0.000214059166, 0.000214059166, 0.258968176, 0.000214059166,
      0.000214059166, 0.000214059166, 0.599666611, 0.000214059166, 0.000214059166, 0},
     {0.000190616771, 0.0678223351, 0.000190616771, 0.000190616771, 0.000190616771, 0.100087597,
      0.000190616771, 0.000190616771, 0.000190616771, 0.830565134, 0.000190616771, 0},
     {0.000101396092, 0.000101396092, 0.0616730324, 0.000101396092, 0.000101396092, 0.000101396092,
      0.143805821, 0.000101396092, 0.000101396092, 0.000101396092, 0.793709977, 0},
     {0.000303994885, 0.000303994885, 0.000303994885, 0.1284746, 0.000303994885, 0.000303994885,
      0.000303994885, 0.141454447, 0.000303994885, 0.000303994885, 0.000303994885, 0.727334999},
     {0.0433312739, 3.22994963e-05, 3.22994963e-05, 3.22994963e-05, 0.117471222, 3.22994963e-05,
      3.22994963e-05, 3.22994963e-05, 0.838939108, 3.22994963e-05, 3.22994963e-05, 0},
     {4.12209605e-05, 0.0329945236, 4.12209605e-05, 4.12209605e-05, 4.12209605e-05, 0.0452292927,
      4.12209605e-05, 4.12209605e-05, 4.12209605e-05, 0.921446416, 4.12209605e-05, 0},
     {5.6318833e-05, 5.6318833e-05, 0.0532155899, 5.6318833e-05, 5.6318833e-05, 5.6318833e-05,
      0.12387295, 5.6318833e-05, 5.6318833e-05, 5.6318833e-05, 0.822460909, 0},
     {0.000171851328, 0.000171851328, 0.000171851328, 0.0769974228, 0.000171851328, 0.000171851328,
      0.000171851328, 0.115058917, 0.000171851328, 0.000171851328, 0.000171851328, 0.806396998},
     {0.0289415541, 3.57312452e-05, 3.57312452e-05, 3.57312452e-05, 0.0863096576, 3.57312452e-05,
      3.57312452e-05, 3.57312452e-05, 0.884462938, 3.57312452e-05, 3.57312452e-05, 0},
     {7.76665226e-05, 0.0487889535, 7.76665226e-05, 7.76665226e-05, 7.76665226e-05, 0.0563094047,
      7.76665226e-05, 7.76665226e-05, 7.76665226e-05, 0.89428031, 7.76665226e-05, 0},
     {0.000224675996, 0.000224675996, 0.11442169, 0.000224675996, 0.000224675996, 0.000224675996,
      0.177555734, 0.000224675996, 0.000224675996, 0.000224675996, 0.706225168, 0},
     {0.000147262818, 0.000147262818, 0.000147262818, 0.0548667631, 0.000147262818, 0.000147262818,
      0.000147262818, 0.0766653324, 0.000147262818, 0.000147262818, 0.000147262818, 0.867142539},
     {0.0251062443, 3.11750569e-05, 3.11750569e-05, 3.11750569e-05, 0.0860636051, 3.11750569e-05,
      3.11750569e-05, 3.11750569e-05, 0.88858075, 3.11750569e-05, 3.11750569e-05, 0},
     {5.64847834e-05, 0.0293039817, 5.64847834e-05, 5.64847834e-05, 5.64847834e-05, 0.0459491268,
      5.64847834e-05, 5.64847834e-05, 5.64847834e-05, 0.924295013, 5.64847834e-05, 0},
     {5.20954226e-05, 5.20954226e-05, 0.0339074593, 5.20954226e-05, 5.20954226e-05, 5.20954226e-05,
      0.103213362, 5.20954226e-05, 5.20954226e-05, 5.20954226e-05, 0.862462415, 0},
     {2.41719133e-05, 2.41719133e-05, 2.41719133e-05, 0.0373123704, 2.41719133e-05, 2.41719133e-05,
      2.41719133e-05, 0.0666923249, 2.41719133e-05, 2.41719133e-05, 2.41719133e-05, 0.895777757}},
    {// stickPmf
     {0.000146474028, 0.0297946523, 0.0803446626, 0.0323153575, 0.000146474028, 0.0181368779,
      0.159977298, 0.0152271765, 0.000146474028, 0.244223284, 0.299539822, 0.120001447},
     {0.0754353435, 0.000106711295, 0.108178173, 0.05166508, 0.0659944961, 0.000106711295,
      0.235061774, 0.0369961862, 0.0285965147, 0.000106711295, 0.241613261, 0.156139038},
     {0.150889801, 0.0424857824, 0.000318026612, 0.147445837, 0.0393052032, 0.0282810617,
      0.000318026612, 0.0807560198, 0.00326633483, 0.214809655, 0.000318026612, 0.291806225},
     {0.140090388, 0.0466498235, 0.0955223418, 0.000146946243, 0.100971439, 0.0170277817,
      0.173755469, 0.000146946243, 0.00152457994, 0.141450628, 0.282713657, 0},
     {0.000119247588, 0.0383163304, 0.11255615, 0.0436285288, 0.000119247588, 0.0192856198,
      0.175758018, 0.020358817, 0.000119247588, 0.197235316, 0.288470104, 0.104033374},
     {0.0468282051, 3.20346125e-05, 0.0825016476, 0.0378754242, 0.0868836623, 3.20346125e-05,
      0.146909879, 0.0191351973, 0.260577625, 3.20346125e-05, 0.213141376, 0.10605088},
     {0.146293011, 0.0331663512, 0.000105430841, 0.0501755794, 0.207187242, 0.0102551561,
      0.000105430841, 0.027888217, 0.283407964, 0.0599407043, 0.000105430841, 0.181369481},
     {0.125236617, 0.0361508505, 0.0591068896, 0.000121340743, 0.165565282, 0.00791956686,
      0.150064084, 0.000121340743, 0.202982991, 0.0283175644, 0.224413472, 0},
     {0.000210016244, 0.0414240647, 0.101760633, 0.0428484744, 0.000210016244, 0.0179312948,
      0.161526364, 0.0189973372, 0.000210016244, 0.240930978, 0.173671163, 0.200279644},
     {0.0495830034, 4.61487199e-05, 0.151014753, 0.067951325, 0.0907268765, 4.61487199e-05,
      0.315461425, 0.0277696796, 0.128849485, 4.61487199e-05, 0.0567220636, 0.111782943},
     {0.11555575, 0.0292950248, 0.00014832042, 0.0456759293, 0.15337244, 0.00761339039,
      0.00014832042, 0.0203902102, 0.249082062, 0.122098355, 0.00014832042, 0.256471877},
     {0.147130835, 0.0327619952, 0.0756043295, 0.000188145734, 0.212444944, 0.0101329647,
      0.114642199, 0.000188145734, 0.294425221, 0.107175814, 0.00530540568, 0},
     {8.22405627e-05, 0.0221237714, 0.037137767, 0.0162121049, 8.22405627e-05, 0.0187609213,
      0.080696449, 0.0269716446, 8.22405627e-05, 0.248406873, 0.205991108, 0.343452639},
     {0.049159614, 5.43832178e-05, 0.0514031038, 0.0203847824, 0.0855463843, 5.43832178e-05,
      0.0980144614, 0.0239136025, 0.168927456, 5.43832178e-05, 0.158422773, 0.344064673},
     {0.137668134, 0.0267105699, 0.000142351166, 0.0441201525, 0.19166368, 0.0114882132,
      0.000142351166, 0.0202277277, 0.273478098, 0.129199413, 0.000142351166, 0.165016958},
     {0.0443720703, 0.011998886, 0.0370775788, 3.30887484e-05, 0.0768170894, 0.00839411198,
      0.110338914, 3.30887484e-05, 0.22919549, 0.235978299, 0.245761384, 0}}};

constexpr double transProbs[nContexts][3][4] = {
    {// AA
     {-7.36174689, 1.61101289, -0.268464948, 0.0143320448},
     {-3.22220388, 0.496986063, -0.181868245, 0.0125555841},
     {2.63396482, -1.49671632, 0.145790372, -0.00463126882}},
    {// AC
     {-4.29124266, 0.474781729, -0.0691005634, 0.0024944116},
     {-5.8759345, 0.769176344, -0.0794793302, 0.0026603066},
     {-2.52278775, -0.0673361054, -0.00420963573, 0.000173101656}},
    {// AG
     {1.44776405, -1.89806596, 0.230489754, -0.00788373432},
     {3.62328388, -3.16155741, 0.380413687, -0.0139731804},
     {-2.4113731, 1.5617699, -0.464222808, 0.0330184761}},
    {// AT
     {-1.52188949, -0.619888164, 0.0344820232, -0.000466695436},
     {7.25871708, -2.89355725, 0.226545422, -0.00503769654},
     {1.80689962, -0.775311988, 0.00287885444, 0.00212630128}},
    {// CA
     {-1.58612186, -0.153996667, -0.019362539, 0.00340134003},
     {-5.44435394, 1.26703171, -0.267105283, 0.0166778421},
     {2.91356449, -1.87150825, 0.226540613, -0.010786308}},
    {// CC
     {2.73738367, -1.52165103, 0.134485954, -0.00364414263},
     {-1.96380465, -0.0304224195, -0.00741383211, 0.000655817377},
     {-5.65275902, 0.8554781, -0.0747563725, 0.00213841659}},
    {// CG
     {6.82408162, -5.00340617, 0.843280331, -0.0477216713},
     {0.922836644, -1.39570517, 0.0615167062, 0.00625866765},
     {-7.49471908, 4.08850937, -0.931670394, 0.06135184}},
    {// CT
     {-2.03471379, -0.183049895, -0.0161448893, 0.00147588228},
     {-0.39397695, -0.497506148, -0.0025279351, 0.00199848358},
     {-3.14511306, 0.866929002, -0.161209014, 0.00724003881}},
    {// GA
     {-4.03613088, 0.951183897, -0.1707077, 0.00951914769},
     {-2.73727929, -0.0846264685, -0.0762702838, 0.00765998521},
     {-0.0115161507, -0.676772845, 0.0186902046, 0.00114061014}},
    {// GC
     {-4.0450919, 0.19441825, -0.0196572463, 0.000545513282},
     {-0.590738944, -0.563607207, 0.0397364279, -0.000835263676},
     {-0.650606525, -0.793978644, 0.0672391429, -0.0019742867}},
    {// GG
     {-7.1216888, 0.704299447, -0.00675682556, -0.00429905236},
     {4.0135083, -2.99745283, 0.330100699, -0.00866415885},
     {-2.95159941, 1.52304256, -0.417790358, 0.0301178899}},
    {// GT
     {-5.71687526, 1.02827154, -0.128832359, 0.00416099843},
     {-6.35542589, 1.25709282, -0.1917998, 0.00889679996},
     {-1.20027791, 0.152425947, -0.0921121841, 0.00518033583}},
    {// TA
     {-1.75233077, -0.0165381386, -0.0154477996, 0.00191469353},
     {-3.79457085, 0.695703623, -0.170024823, 0.0125082712},
     {2.14044448, -1.41498734, 0.116355694, -0.00281465688}},
    {// TC
     {-1.89375385, -0.103804215, -0.00397067386, 0.000435442231},
     {-4.05701629, 0.431318996, -0.0460802272, 0.0014783096},
     {-4.07119813, 0.0340792867, 0.00636890829, -0.000665968456}},
    {// TG
     {7.14617037, -5.07757031, 0.837708073, -0.0452627753},
     {1.69236944, -2.28936276, 0.274337077, -0.00888968456},
     {-4.60392728, 2.15828049, -0.515604387, 0.0324594302}},
    {// TT
     {-2.82895825, -0.159891687, 0.0556441946, -0.00321555705},
     {2.92995068, -1.83116699, 0.202125712, -0.00725218663},
     {0.981168958, -0.682236621, 0.0147684255, 0.00146373056}}};

inline double CalculateExpectedLLForEmission(const size_t move, const uint8_t row,
                                             const size_t moment)
{
    double expectedLL = 0;
    for (size_t i = 0; i < nOutcomes; i++) {
        double curProb = emissionPmf[move][row][i];
        double lgCurProb = std::log(curProb);
        if (moment == static_cast<uint8_t>(MomentType::FIRST))
            expectedLL += curProb * lgCurProb;
        else if (moment == static_cast<uint8_t>(MomentType::SECOND))
            expectedLL += curProb * (lgCurProb * lgCurProb);
    }
    return expectedLL;
}

S_P2C2_Model::S_P2C2_Model(const SNR& snr) : snr_(snr)
{
    for (size_t ctx = 0; ctx < nContexts; ++ctx) {
        const uint8_t bp = ctx & 3;  // (equivalent to % 4)
        const double snr1 = clip(snr_[bp], snrRanges[bp]), snr2 = snr1 * snr1, snr3 = snr2 * snr1;
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

std::unique_ptr<AbstractRecursor> S_P2C2_Model::CreateRecursor(
    std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr, double scoreDiff) const
{
    return std::unique_ptr<AbstractRecursor>(
        new S_P2C2_Recursor(std::forward<std::unique_ptr<AbstractTemplate>>(tpl), mr, scoreDiff));
}

std::vector<TemplatePosition> S_P2C2_Model::Populate(const std::string& tpl) const
{
    std::vector<TemplatePosition> result;

    if (tpl.empty()) return result;

    result.reserve(tpl.size());

    // calculate transition probabilities
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

double S_P2C2_Model::ExpectedLLForEmission(const MoveType move, const uint8_t prev,
                                           const uint8_t curr, const MomentType moment) const
{
    const size_t row = (prev << 2) | curr;
    return cachedEmissionExpectations_[row][static_cast<uint8_t>(move)]
                                      [static_cast<uint8_t>(moment)];
}

S_P2C2_Recursor::S_P2C2_Recursor(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr,
                                 double scoreDiff)
    : Recursor<S_P2C2_Recursor>(std::forward<std::unique_ptr<AbstractTemplate>>(tpl), mr, scoreDiff)
{
}

std::vector<uint8_t> S_P2C2_Recursor::EncodeRead(const MappedRead& read)
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

double S_P2C2_Recursor::EmissionPr(MoveType move, uint8_t emission, uint8_t prev, uint8_t curr)
{
    assert(move != MoveType::DELETION);
    const auto row = (prev << 2) | curr;
    return emissionPmf[static_cast<uint8_t>(move)][row][emission] * kCounterWeight;
}

double S_P2C2_Recursor::UndoCounterWeights(const size_t nEmissions) const
{
    return -std::log(kCounterWeight) * nEmissions;
}
}  // namespace anonymous
}  // namespace Consensus
}  // namespace PacBio
