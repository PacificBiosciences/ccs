
#include "Mutations.h"

std::vector<PacBio::Consensus::Mutation> Mutations(const std::string& tpl, const size_t start,
                                                   const size_t end)
{
    using namespace PacBio::Consensus;

    constexpr auto bases = "ACGT";

    std::vector<Mutation> result;

    for (size_t i = start; i < end; ++i) {
        for (size_t j = 0; j < 4; ++j)
            result.push_back(Mutation(MutationType::INSERTION, i, bases[j]));

        result.push_back(Mutation(MutationType::DELETION, i));

        for (size_t j = 0; j < 4; ++j)
            if (bases[j] != tpl[i])
                result.push_back(Mutation(MutationType::SUBSTITUTION, i, bases[j]));
    }

    for (size_t j = 0; j < 4; ++j)
        result.push_back(Mutation(MutationType::INSERTION, tpl.length(), bases[j]));

    return result;
}

std::vector<PacBio::Consensus::Mutation> Mutations(const std::string& tpl)
{
    return Mutations(tpl, 0, tpl.length());
}
