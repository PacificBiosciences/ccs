// Author: Lance Hepler

#include <algorithm>
#include <cmath>
#include <functional>
#include <iterator>
#include <limits>
#include <list>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/optional.hpp>

#include <pbcopper/logging/Logging.h>

#include <pacbio/consensus/Integrator.h>
#include <pacbio/consensus/Polish.h>
#include <pacbio/exception/InvalidEvaluatorException.h>

#include "MutationTracker.h"

using std::list;
using std::pair;
using std::set;
using std::string;
using std::vector;

using std::make_pair;
using std::tie;

namespace PacBio {
namespace Consensus {

PolishConfig::PolishConfig(const size_t iterations, const size_t separation,
                           const size_t neighborhood, const bool diploid)
    : MaximumIterations(iterations)
    , MutationSeparation(separation)
    , MutationNeighborhood(neighborhood)
    , Diploid(diploid)
{
}

RepeatConfig::RepeatConfig(const size_t repeatSize, const size_t elementCount,
                           const size_t iterations)
    : MaximumRepeatSize{repeatSize}
    , MinimumElementCount{elementCount}
    , MaximumIterations{iterations}
{
}

void Mutations(vector<Mutation>* muts, const Integrator& ai, const size_t start, const size_t end,
               const bool diploid = false)
{
    const std::vector<char> bases{
        diploid ? std::vector<char>{'A', 'C', 'G', 'T', 'Y', 'R', 'W', 'S', 'K', 'M'}
                : std::vector<char>{'A', 'C', 'G', 'T'}};

    std::function<bool(const char&, const char&)> containedWithin;
    if (diploid) {
        // in diploid mode, we want to generate candidates
        // that are unequal to the current *char*, i.e.,
        // say we have a 'Y' (='C'+'T'), we still want to
        // generate a 'C' and 'T', we just don't want a 'Y'.
        containedWithin = std::equal_to<char>{};
    } else {
        // in haploid mode, we want to avoid all *subsets*
        // of pure bases, that is, if curr is 'Y', then we
        // neither want a 'C' nor a 'T'.
        containedWithin = Data::detail::ambiguousBaseContainsPureBase;
    }

    if (start == end) return;

    char last = (start > 0) ? ai[start - 1] : '\0';

    for (size_t i = start; i < end; ++i) {
        const char curr = ai[i];

        // insertions come before deletion/substitutions at site i, their End()
        // is i < i + 1
        for (const char j : bases) {
            // if it's not a homopolymer insertion, or it's the first base of
            // one..
            if (!containedWithin(last, j)) muts->emplace_back(Mutation::Insertion(i, j));
        }

        // if we're the first in the homopolymer, we can delete
        if (curr != last) muts->emplace_back(Mutation::Deletion(i, 1));

        for (const char j : bases) {
            if (!containedWithin(curr, j)) muts->emplace_back(Mutation::Substitution(i, j));
        }

        last = curr;
    }

    // if we are at the end, make sure we're not performing a terminal
    // homopolymer insertion
    for (const char j : bases)
        if (!containedWithin(last, j)) muts->emplace_back(Mutation::Insertion(end, j));
}

vector<Mutation> Mutations(const Integrator& ai, const size_t start, const size_t end,
                           const bool diploid = false)
{
    vector<Mutation> muts;
    Mutations(&muts, ai, start, end, diploid);
    return muts;
}

vector<Mutation> Mutations(const Integrator& ai, const bool diploid)
{
    return Mutations(ai, 0, ai.TemplateLength(), diploid);
}

void RepeatMutations(vector<Mutation>* muts, const Integrator& ai, const RepeatConfig& cfg,
                     const size_t start, const size_t end)
{
    if (cfg.MaximumRepeatSize < 2 || cfg.MinimumElementCount <= 0) return;

    const string tpl(ai);

    for (size_t repeatSize = 2; repeatSize <= cfg.MaximumRepeatSize; ++repeatSize) {
        for (size_t i = start; i + repeatSize <= end;) {
            size_t nElem = 1;

            for (size_t j = i + repeatSize; j + repeatSize <= end; j += repeatSize) {
                if (tpl.compare(j, repeatSize, tpl, i, repeatSize) == 0)
                    ++nElem;
                else
                    break;
            }

            if (nElem >= cfg.MinimumElementCount) {
                muts->emplace_back(Mutation::Insertion(i, tpl.substr(i, repeatSize)));
                muts->emplace_back(Mutation::Deletion(i, repeatSize));
            }

            if (nElem > 1)
                i += repeatSize * (nElem - 1) + 1;
            else
                ++i;
        }
    }

    sort(muts->begin(), muts->end(), Mutation::SiteComparer);
}

vector<Mutation> RepeatMutations(const Integrator& ai, const RepeatConfig& cfg, const size_t start,
                                 const size_t end)
{
    vector<Mutation> muts;
    RepeatMutations(&muts, ai, cfg, start, end);
    return muts;
}

vector<Mutation> RepeatMutations(const Integrator& ai, const RepeatConfig& cfg)
{
    return RepeatMutations(ai, cfg, 0, ai.TemplateLength());
}

vector<Mutation> BestMutations(list<ScoredMutation>* scoredMuts, const size_t separation)
{
    vector<Mutation> result;

    // TODO handle 0-separation correctly
    if (separation == 0) throw std::invalid_argument("nonzero separation required");

    while (!scoredMuts->empty()) {
        const auto& mut =
            *max_element(scoredMuts->begin(), scoredMuts->end(), ScoredMutation::ScoreComparer);

        result.emplace_back(mut);

        const size_t start = (separation < mut.Start()) ? mut.Start() - separation : 0;
        const size_t end = mut.End() + separation;

        scoredMuts->remove_if(
            [start, end](const ScoredMutation& m) { return start <= m.End() && m.Start() < end; });
    }

    return result;
}

vector<Mutation> NearbyMutations(vector<Mutation>* applied, vector<Mutation>* centers,
                                 const Integrator& ai, const size_t neighborhood,
                                 const bool diploid = false)
{
    const size_t len = ai.TemplateLength();
    const auto clamp = [len](const int i) { return std::max(0, std::min<int>(len, i)); };

    vector<Mutation> result;

    if (centers->empty()) return result;

    sort(applied->begin(), applied->end(), Mutation::SiteComparer);
    sort(centers->begin(), centers->end(), Mutation::SiteComparer);

    const auto mutRange = [clamp, neighborhood](const Mutation& mut, const int diff) {
        const int start = diff + mut.Start() - neighborhood;
        const int end = diff + mut.End() + neighborhood;
        return pair<size_t, size_t>(clamp(start), clamp(end));
    };

    // find the ranges
    auto ait = applied->cbegin();
    auto cit = centers->cbegin();
    int lengthDiff = 0;

    for (; ait != applied->cend() && ait->End() <= cit->Start(); ++ait)
        lengthDiff += ait->LengthDiff();

    vector<pair<size_t, size_t>> ranges = {mutRange(*cit, lengthDiff)};
    size_t currEnd = ranges.back().second;

    // increment to the next centerpoint and continue
    for (++cit; cit != centers->cend(); ++cit) {
        size_t nextStart, nextEnd;

        for (; ait != applied->cend() && ait->End() <= cit->Start(); ++ait)
            lengthDiff += ait->LengthDiff();

        tie(nextStart, nextEnd) = mutRange(*cit, lengthDiff);

        // if the next range touches the last one, just extend the last one
        if (nextStart <= currEnd)
            ranges.back().second = nextEnd;
        else {
            ranges.emplace_back(make_pair(nextStart, nextEnd));
            currEnd = nextEnd;
        }
    }

    for (const auto& range : ranges)
        Mutations(&result, ai, range.first, range.second, diploid);

    return result;
}

// The significance level for the likelihood-ratio test of
// rejecting the null of having a purely haploid site.
// We use 0.5%, in order to make strong claims for our discoveries
// Reference:
//   https://www.nature.com/articles/s41562-017-0189-z
static constexpr const double significanceLevel = 0.005;

PolishResult Polish(Integrator* ai, const PolishConfig& cfg)
{
    vector<Mutation> muts = Mutations(*ai, cfg.Diploid);
    std::hash<string> hashFn;
    size_t oldTpl = hashFn(*ai);
    set<size_t> history = {oldTpl};

    PolishResult result;
    // keep track of the changes to the original template over many rounds
    MutationTracker mutTracker{static_cast<std::string>(*ai)};

    // In Haploid mode, the LL just needs to improve, i.e.,
    //
    //   newLL > currentLL
    //
    // or equivalently
    //
    //   newLL - currentLL > 0
    //
    // In Diploid mode however, we perform an implicit likelihood
    // ratio test, where
    //
    //   H_0: current haploid base
    //   H_A: prospective diploid base
    //
    // We make an extremely conservative assumption that H_A has
    // 3 degrees of freedom more than H_0 (which is not true, but
    // only makes the test more conservative, i.e., we trade a
    // [negligible] amount of sensitivity for more specificity).
    // To calculate this we have to calculate the ChiSq quantile
    // function for (1 - significanceLevel).
    const double minImprovementThreshold{
        cfg.Diploid
            ? boost::math::quantile(complement(boost::math::chi_squared{3}, significanceLevel))
            : 0};

    for (size_t i = 0; i < cfg.MaximumIterations; ++i) {
        // find the best mutations given our parameters
        {
            list<ScoredMutation> scoredMuts;
            int mutationsTested = 0;
            bool hasNewInvalidEvaluator;

            // Compute new sets of possible mutations until no Evaluators are
            // being invalidated.
            do {
                // Compute the LL only with the active Evaluators
                const double LL = ai->LL();

                hasNewInvalidEvaluator = false;
                try {
                    // Get set of possible mutations
                    for (const auto& mut : muts) {
                        ++mutationsTested;
                        const double ll = ai->LL(mut);
                        if (ll - LL > (mut.IsDeletion() ? 0 : minImprovementThreshold))
                            scoredMuts.emplace_back(mut.WithScore(ll));
                    }
                } catch (const Exception::InvalidEvaluatorException& e) {
                    // If an Evaluator exception occured,
                    // retry without problematic Evaluator
                    PBLOG_INFO << e.what();
                    hasNewInvalidEvaluator = true;
                    scoredMuts.clear();
                    mutationsTested = 0;
                }
            } while (hasNewInvalidEvaluator);

            result.mutationsTested += mutationsTested;

            // take best mutations in separation window, apply them
            muts = BestMutations(&scoredMuts, cfg.MutationSeparation);
        }

        // convergence!!
        if (muts.empty()) {
            result.hasConverged = true;

            if (cfg.Diploid) {
                result.diploidSites = mutTracker.MappingToOriginalTpl();
            }

            return result;
        }

        const size_t newTpl = hashFn(ApplyMutations(*ai, &muts));

        if (cfg.Diploid) {
            mutTracker.AddSortedMutations(muts);
        }

        const auto diagnostics = [&result](Integrator* ai) {
            result.maxAlphaPopulated.emplace_back(ai->MaxAlphaPopulated());
            result.maxBetaPopulated.emplace_back(ai->MaxBetaPopulated());
            result.maxNumFlipFlops.emplace_back(ai->MaxNumFlipFlops());
        };

        if (history.find(newTpl) != history.end()) {
            /* Cyclic behavior guard - Dave A. found some edge cases where the
             template was mutating back to an earlier version. This is a bad
             and should be rare.  He found that by applying the single best
             mutation you could avoid the loop. (That is if adding Muts X + Y
             made removing muts X + Y beneficial, then you can break that
             inifinite loop by just applying X or Y, as presumably this removes
             the interaction between them that leads to the cycling behavior.
             This step is just a heuristic work around that was found. */
            ai->ApplyMutation(muts.front());
            oldTpl = hashFn(*ai);
            ++result.mutationsApplied;

            diagnostics(ai);

            // get the mutations for the next round
            vector<Mutation> applied = {muts.front()};
            muts = NearbyMutations(&applied, &muts, *ai, cfg.MutationNeighborhood, cfg.Diploid);
        } else {
            ai->ApplyMutations(&muts);
            oldTpl = newTpl;
            result.mutationsApplied += muts.size();

            diagnostics(ai);

            // get the mutations for the next round
            muts = NearbyMutations(&muts, &muts, *ai, cfg.MutationNeighborhood, cfg.Diploid);
        }

        // keep track of which templates we've seen
        history.insert(oldTpl);
    }

    return result;
}

PolishResult PolishRepeats(Integrator* const ai, const RepeatConfig& cfg)
{
    PolishResult result;

    const auto diagnostics = [&result](Integrator* ai) {
        result.maxAlphaPopulated.emplace_back(ai->MaxAlphaPopulated());
        result.maxBetaPopulated.emplace_back(ai->MaxBetaPopulated());
        result.maxNumFlipFlops.emplace_back(ai->MaxNumFlipFlops());
    };

    for (size_t i = 0; i < cfg.MaximumIterations; ++i) {
        const vector<Mutation> muts = RepeatMutations(*ai, cfg);
        boost::optional<ScoredMutation> bestMut = boost::none;
        size_t mutationsTested = 0;
        bool hasNewInvalidEvaluator = false;

        // if an Evaluator exception occurs, restart
        do {
            const double LL = ai->LL();
            hasNewInvalidEvaluator = false;
            try {
                for (const auto& mut : muts) {
                    const double ll = ai->LL(mut);
                    if (ll > LL && (!bestMut || bestMut->Score < ll)) bestMut = mut.WithScore(ll);
                }
            } catch (const Exception::InvalidEvaluatorException& e) {
                PBLOG_INFO << e.what();
                hasNewInvalidEvaluator = true;
                bestMut = boost::none;
                mutationsTested = 0;
            }
        } while (hasNewInvalidEvaluator);

        result.mutationsTested += mutationsTested;

        if (!bestMut) {
            result.hasConverged = true;
            break;
        }

        std::vector<Mutation> mut = {Mutation(*bestMut)};
        ai->ApplyMutations(&mut);
        ++result.mutationsApplied;
        diagnostics(ai);
    }

    return result;
}

namespace {  // anonymous

int ProbabilityToQV(double probability)
{
    if (probability < 0.0 || probability > 1.0)
        throw std::invalid_argument("invalid value: probability not in [0,1]");
    else if (probability == 0.0)
        probability = std::numeric_limits<double>::min();

    return static_cast<int>(round(-10.0 * log10(probability)));
}

inline int ScoreSumToQV(const double scoreSum)
{
    return ProbabilityToQV(1.0 - 1.0 / (1.0 + scoreSum));
}

}  // anonymous namespace

vector<int> ConsensusQualities(Integrator& ai)
{
    vector<int> quals;
    quals.reserve(ai.TemplateLength());
    const double LL = ai.LL();
    for (size_t i = 0; i < ai.TemplateLength(); ++i) {
        double scoreSum = 0.0;
        for (const auto& m : Mutations(ai, i, i + 1)) {
            // skip mutations that start beyond the current site (e.g. trailing insertions)
            if (m.Start() > i) continue;
            // TODO (lhepler): this is dumb, but untestable mutations,
            //   aka insertions at ends, cause all sorts of weird issues
            // See also: Polish::ConsensusQVs(ai)
            double score;
            try {
                score = ai.LL(m) - LL;
            } catch (const Exception::InvalidEvaluatorException& e) {
                // If an Evaluator exception occured, report and skip!
                // We need to handle this!
                std::string error = "In Polish::ConsensusQualities(ai): ";
                error += e.what();
                PBLOG_ERROR << error;
                continue;
            }
            assert(score <= 0.0);

            if (score < 0) scoreSum += exp(score);
        }
        quals.emplace_back(ScoreSumToQV(scoreSum));
    }
    return quals;
}

QualityValues ConsensusQVs(Integrator& ai)
{
    const size_t len = ai.TemplateLength();
    vector<int> quals, delQVs, insQVs, subQVs;
    quals.reserve(len);
    delQVs.reserve(len);
    insQVs.reserve(len);
    subQVs.reserve(len);
    const double LL = ai.LL();
    for (size_t i = 0; i < len; ++i) {
        double qualScoreSum = 0.0, delScoreSum = 0.0, insScoreSum = 0.0, subScoreSum = 0.0;
        for (const auto& m : Mutations(ai, i, i + 1)) {
            // skip mutations that start beyond the current site (e.g. trailing insertions)
            if (m.Start() > i) continue;

            // TODO (lhepler): this is dumb, but untestable mutations,
            //   aka insertions at ends, cause all sorts of weird issues
            // See also: Polish::ConsensusQualities(ai)
            double score;
            try {
                score = ai.LL(m) - LL;
            } catch (const Exception::InvalidEvaluatorException& e) {
                // If an Evaluator exception occured, report and skip!
                // We need to handle this!
                std::string error = "In Polish::ConsensusQVs(ai): ";
                error += e.what();
                PBLOG_ERROR << error;
                continue;
            }

            // this really should never happen
            if (score >= 0.0) continue;
            const double expScore = exp(score);
            qualScoreSum += expScore;
            if (m.IsDeletion())
                delScoreSum += expScore;
            else if (m.Start() == m.End())
                insScoreSum += expScore;
            else
                subScoreSum += expScore;
        }
        quals.emplace_back(ScoreSumToQV(qualScoreSum));
        delQVs.emplace_back(ScoreSumToQV(delScoreSum));
        insQVs.emplace_back(ScoreSumToQV(insScoreSum));
        subQVs.emplace_back(ScoreSumToQV(subScoreSum));
    }
    // TODO(lhepler): discuss InsQV being len + 1 to capture trailing insertion
    return QualityValues{std::move(quals), std::move(delQVs), std::move(insQVs), std::move(subQVs)};
}

}  // namespace Consensus
}  // namespace PacBio
