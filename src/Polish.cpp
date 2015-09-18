
#include <algorithm>
#include <functional>
#include <iterator>
#include <list>
#include <set>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include <pacbio/consensus/Integrator.h>
#include <pacbio/consensus/Polish.h>

namespace PacBio {
namespace Consensus {

PolishConfig::PolishConfig(const size_t iterations, const size_t separation,
                           const size_t neighborhood)
    : MaximumIterations(iterations)
    , MutationSeparation(separation)
    , MutationNeighborhood(neighborhood)
{
}

void Mutations(std::vector<Mutation>* muts, const AbstractIntegrator& ai, const size_t start,
               const size_t end)
{
    constexpr auto bases = "ACGT";

    if (start == end) return;

    char last = (start > 0) ? ai[start - 1] : '\0';

    for (size_t i = start; i < end; ++i) {
        const char curr = ai[i];

        // insertions come before deletion/substitutions at site i, their End()
        // is i < i + 1
        for (size_t j = 0; j < 4; ++j) {
            // if it's not a homopolymer insertion, or it's the first base of
            // one..
            if (bases[j] != last)
                muts->emplace_back(Mutation(MutationType::INSERTION, i, bases[j]));
        }

        // if we're the first in the homopolymer, we can delete
        if (curr != last) muts->emplace_back(Mutation(MutationType::DELETION, i));

        for (size_t j = 0; j < 4; ++j) {
            if (bases[j] != curr)
                muts->emplace_back(Mutation(MutationType::SUBSTITUTION, i, bases[j]));
        }

        last = curr;
    }

    // if we're not at the absolute end, use prior algorithm
    if (end < ai.Length()) {
        const char curr = ai[end];

        for (size_t j = 0; j < 4; ++j)
            if (bases[j] != last)
                muts->emplace_back(Mutation(MutationType::INSERTION, end, bases[j]));
    }
    // if we are at the end, make sure we're not performing a terminal
    // homopolymer insertion
    else {
        for (size_t j = 0; j < 4; ++j)
            if (bases[j] != last)
                muts->emplace_back(Mutation(MutationType::INSERTION, end, bases[j]));
    }
}

std::vector<Mutation> Mutations(const AbstractIntegrator& ai, const size_t start, const size_t end)
{
    std::vector<Mutation> muts;
    Mutations(&muts, ai, start, end);
    return muts;
}

std::vector<Mutation> Mutations(const AbstractIntegrator& ai)
{
    return Mutations(ai, 0, ai.Length());
}

std::vector<Mutation> BestMutations(std::list<ScoredMutation>* scoredMuts, const size_t separation)
{
    std::vector<Mutation> muts;

    if (separation == 0) {
        std::copy(scoredMuts->begin(), scoredMuts->end(), std::back_inserter(muts));
        return muts;
    }

    while (!scoredMuts->empty()) {
        const auto& mut = *std::max_element(scoredMuts->begin(), scoredMuts->end(),
                                            ScoredMutation::ScoreComparer);

        muts.emplace_back(mut);

        const size_t start = (separation < mut.Start()) ? mut.Start() - separation : 0;
        const size_t end = mut.End() + separation;

        scoredMuts->remove_if(
            [start, end](const ScoredMutation& m) { return start <= m.End() && m.Start() < end; });
    }

    return muts;
}

std::vector<Mutation> NearbyMutations(const AbstractIntegrator& ai,
                                      const std::vector<Mutation>& centers,
                                      const size_t neighborhood)
{
    const size_t len = ai.Length();

    std::vector<Mutation> muts;

    if (centers.empty()) return muts;

    auto mutRange = [len, neighborhood](const Mutation& mut) {
        const size_t start = (neighborhood > mut.Start()) ? 0 : mut.Start() - neighborhood;
        const size_t end = std::min(len, mut.End() + neighborhood);
        return std::make_tuple(start, end);
    };

    // find the ranges
    std::vector<Mutation>::const_iterator it = centers.begin();
    std::vector<std::tuple<size_t, size_t>> ranges = {mutRange(*it)};
    std::tuple<size_t, size_t>& lastRange = ranges.back();

    // increment to the next position and continue
    for (++it; it != centers.end(); ++it) {
        const auto nextRange = mutRange(*it);
        // if the new range touches the last one, just extend the last one
        if (std::get<0>(nextRange) <= std::get<1>(lastRange))
            std::get<1>(lastRange) = std::get<1>(nextRange);
        else {
            ranges.push_back(nextRange);
            lastRange = ranges.back();
        }
    }

    for (const auto& range : ranges)
        Mutations(&muts, ai, std::get<0>(range), std::get<1>(range));

    return muts;
}

bool Polish(AbstractIntegrator* ai, const PolishConfig& cfg)
{
    std::vector<Mutation> muts = Mutations(*ai);
    std::hash<std::string> hashFn;
    size_t oldTpl = hashFn(*ai);
    std::set<size_t> history = {oldTpl};

    for (size_t i = 0; i < cfg.MaximumIterations; ++i) {
        boost::optional<Mutation> bestMut;

        // find the best mutations given our parameters
        {
            const double LL = ai->LL();
            std::list<ScoredMutation> scoredMuts;
            double bestll = LL;

            for (const auto& mut : muts) {
                const double ll = ai->LL(mut);
                if (ll > LL) scoredMuts.emplace_back(mut.WithScore(ll));
                if (ll > bestll) bestMut = mut;
            }

            // take best mutations in separation window, apply them
            muts = BestMutations(&scoredMuts, cfg.MutationSeparation);
        }

        const size_t newTpl = hashFn(ApplyMutations(*ai, &muts));

        // convergence!!
        if (newTpl == oldTpl) return true;

        if (history.find(newTpl) != history.end()) {
            // cyclic behavior detected! apply just the single best mutation
            if (!bestMut)
                throw std::runtime_error(
                    "entered cycle detection without any mutations to test...");

            ai->ApplyMutation(*bestMut);
            oldTpl = hashFn(*ai);
        } else {
            ai->ApplyMutations(&muts);
            oldTpl = newTpl;
        }

        // keep track of which templates we've seen
        history.insert(oldTpl);

        // get the mutations for the next round
        muts = NearbyMutations(*ai, muts, cfg.MutationNeighborhood);
    }

    return false;
}

}  // namespace Consensus
}  // namespace PacBio
