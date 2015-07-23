
#include <cmath>
#include <utility>

#include <ConsensusCore/Arrow/ContextParameters.hpp>
#include <ConsensusCore/Arrow/TemplateParameterPair.hpp>

namespace ConsensusCore {
namespace Arrow {
namespace {

inline
std::pair<double, double> ExpectedContextLL(const TransitionParameters& params, double eps)
{
    const double p_m = params.Match,     l_m = std::log(p_m),  l2_m = l_m * l_m;
    const double p_d = params.Deletion,  l_d = std::log(p_d),  l2_d = l_d * l_d;
    const double p_b = params.Branch,    l_b = std::log(p_b),  l2_b = l_b * l_b;
    const double p_s = params.Stick,     l_s = std::log(p_s),  l2_s = l_s * l_s;

    const double lgThird = -std::log(3.0);
    const double E_M = (1.0 - eps) * 0.0 + eps * lgThird,  E2_M = E_M * E_M;
    const double E_D = 0.0,                                E2_D = E_D * E_D;
    const double E_B = 0.0,                                E2_B = E_B * E_B;
    const double E_S = lgThird,                            E2_S = E_S * E_S;

    auto ENN = [=] (const double l_m, const double l_d, const double l_b, const double l_s,
                    const double E_M, const double E_D, const double E_B, const double E_S)
    {
        const double E_MD = (l_m + E_M) * p_m / (p_m + p_d) + (l_d + E_D) * p_d / (p_m + p_d);
        const double E_I  = (l_b + E_B) * p_b / (p_b + p_s) + (l_s + E_S) * p_s / (p_b + p_s);
        const double E_BS = E_I * (p_s + p_b) / (p_m + p_d); 
        return E_MD + E_BS;
    };

    const double mean = ENN(l_m, l_d, l_b, l_s, E_M, E_D, E_B, E_S);
    const double var  = ENN(l2_m, l2_d, l2_b, l2_s, E2_M, E2_D, E2_B, E2_S) - mean * mean;

    return std::make_pair(mean, var);
}

}

inline
std::vector<std::pair<double, double>> PerBaseMeanAndVariance(const TemplateParameterPair& tpl,
                                                              double eps)
{
    std::vector<std::pair<double, double>> meanVars;

    for (size_t i = 0; i < tpl.Length(); ++i)
    {
        meanVars.emplace_back(ExpectedContextLL(tpl.GetTemplatePosition(i).second, eps));
    }

    return meanVars;
}

} // namespace Arrow
} // namespace ConsensusCore
