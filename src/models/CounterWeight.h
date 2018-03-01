// Author: Lance Hepler

#pragma once

#include <cmath>
#include <functional>

#include <pacbio/consensus/ModelConfig.h>

namespace PacBio {
namespace Consensus {
namespace {

inline double CounterWeight(std::function<double(size_t, MoveType)>&& lgPrTransition,
                            std::function<double(size_t, MoveType)>&& lgPrEmission,
                            const size_t nContexts)
{
    double meanPrEm = 0.0;
    for (size_t ctx = 0; ctx < nContexts; ++ctx) {
        const double pr_M = lgPrTransition(ctx, MoveType::MATCH);
        const double pr_B = lgPrTransition(ctx, MoveType::BRANCH);
        const double pr_S = lgPrTransition(ctx, MoveType::STICK);
        const double pr_D = lgPrTransition(ctx, MoveType::DELETION);

        const double lgPr_EM = lgPrEmission(ctx, MoveType::MATCH);
        const double lgPr_EB = lgPrEmission(ctx, MoveType::BRANCH);
        const double lgPr_ES = lgPrEmission(ctx, MoveType::STICK);
        const double lgPr_ED = 0.0;  // nothing to emit

        const double E_MD = lgPr_EM * pr_M / (pr_M + pr_D) + lgPr_ED * pr_D / (pr_M + pr_D);
        const double E_BS = lgPr_EB * pr_B / (pr_B + pr_S) + lgPr_ES * pr_S / (pr_B + pr_S);
        const double E_I = E_BS * (pr_B + pr_S) / (pr_M + pr_D);

        meanPrEm += std::exp(E_I + E_MD);
    }
    meanPrEm /= nContexts;
    return 1.0 / meanPrEm;
}
}
}
}
