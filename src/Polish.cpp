
#include <algorithm>
#include <cmath>
#include <functional>
#include <iterator>
#include <list>
#include <limits>
#include <set>
#include <string>
#include <tuple>
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

std::vector<Mutation> NearbyMutations(std::vector<Mutation>* applied,
                                      std::vector<Mutation>* centers, const AbstractIntegrator& ai,
                                      const size_t neighborhood)
{
    const size_t len = ai.Length();
    const auto clamp = [len](const int i) { return std::max(0, std::min<int>(len, i)); };

    std::vector<Mutation> muts;

    if (centers->empty()) return muts;

    std::sort(applied->begin(), applied->end(), Mutation::SiteComparer);
    std::sort(centers->begin(), centers->end(), Mutation::SiteComparer);

    const auto mutRange = [clamp, neighborhood](const Mutation& mut, const int diff) {
        const int start = diff + mut.Start() - neighborhood;
        const int end = diff + mut.End() + neighborhood;
        return std::tuple<size_t, size_t>(clamp(start), clamp(end));
    };

    // find the ranges
    auto ait = applied->cbegin();
    auto cit = centers->cbegin();
    int lengthDiff = 0;

    for (; ait != applied->cend() && ait->End() <= cit->Start(); ++ait)
        lengthDiff += ait->LengthDiff();

    std::vector<std::tuple<size_t, size_t>> ranges = {mutRange(*cit, lengthDiff)};
    size_t currEnd = std::get<1>(ranges.back());

    // increment to the next position and continue
    for (++cit; cit != centers->cend(); ++cit) {
        size_t nextStart, nextEnd;

        for (; ait != applied->cend() && ait->End() <= cit->Start(); ++ait)
            lengthDiff += ait->LengthDiff();

        std::tie(nextStart, nextEnd) = mutRange(*cit, lengthDiff);

        // if the next range touches the last one, just extend the last one
        if (nextStart <= currEnd)
            std::get<1>(ranges.back()) = nextEnd;
        else {
            ranges.emplace_back(std::make_tuple(nextStart, nextEnd));
            currEnd = nextEnd;
        }
    }

    for (const auto& range : ranges)
        Mutations(&muts, ai, std::get<0>(range), std::get<1>(range));

    return muts;
}

std::tuple<bool, size_t, size_t> Polish(AbstractIntegrator* ai, const PolishConfig& cfg)
{
    std::vector<Mutation> muts = Mutations(*ai);
    std::hash<std::string> hashFn;
    size_t oldTpl = hashFn(*ai);
    std::set<size_t> history = {oldTpl};
    size_t nTested = 0;
    size_t nApplied = 0;

    for (size_t i = 0; i < cfg.MaximumIterations; ++i) {
        // find the best mutations given our parameters
        {
            const double LL = ai->LL();
            std::list<ScoredMutation> scoredMuts;

            for (const auto& mut : muts) {
                const double ll = ai->LL(mut);
                if (ll > LL) scoredMuts.emplace_back(mut.WithScore(ll));
                ++nTested;
            }

            // take best mutations in separation window, apply them
            muts = BestMutations(&scoredMuts, cfg.MutationSeparation);
        }

        const size_t newTpl = hashFn(ApplyMutations(*ai, &muts));

        // convergence!!
        if (newTpl == oldTpl) return std::make_tuple(true, nTested, nApplied);

        if (history.find(newTpl) != history.end()) {
            // cyclic behavior detected! apply just the single best mutation
            if (muts.empty())
                throw std::runtime_error(
                    "entered cycle detection without any mutations to test...");

            ai->ApplyMutation(muts.front());
            oldTpl = hashFn(*ai);
            ++nApplied;

            // get the mutations for the next round
            std::vector<Mutation> applied(muts.begin(), muts.begin() + 1);
            muts = NearbyMutations(&applied, &muts, *ai, cfg.MutationNeighborhood);
        } else {
            ai->ApplyMutations(&muts);
            oldTpl = newTpl;
            nApplied += muts.size();

            // get the mutations for the next round
            muts = NearbyMutations(&muts, &muts, *ai, cfg.MutationNeighborhood);
        }

        // keep track of which templates we've seen
        history.insert(oldTpl);
    }

    return std::make_tuple(false, nTested, nApplied);
}

int ProbabilityToQV(double probability)
{
    if (probability < 0.0 || probability > 1.0)
        throw std::invalid_argument("invalid value: probability not in [0,1]");
    else if (probability == 0.0)
        probability = std::numeric_limits<double>::min();

    return static_cast<int>(round(-10.0 * log10(probability)));
}

std::vector<int> ConsensusQVs(AbstractIntegrator& ai)
{
    std::vector<int> QVs;
    const double LL = ai.LL();
    for (size_t i = 0; i < ai.Length(); ++i) {
        double scoreSum = 0.0;
        for (const auto& m : Mutations(ai, i, i + 1)) {
            // skip mutations that start beyond the current site (e.g. trailing insertions)
            if (m.Start() > i) continue;
            // TODO (lhepler): this is dumb, but untestable mutations,
            //   aka insertions at ends, cause all sorts of weird issues
            double score = ai.LL(m) - LL;
            if (score < 0.0) scoreSum += exp(score);
        }
        QVs.push_back(ProbabilityToQV(1.0 - 1.0 / (1.0 + scoreSum)));
    }
    return QVs;
}

}  // namespace Consensus
}  // namespace PacBio
