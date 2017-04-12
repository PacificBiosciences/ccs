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

#include <boost/optional.hpp>

#include <pbcopper/logging/Logging.h>

#include <pacbio/consensus/AbstractIntegrator.h>
#include <pacbio/consensus/Polish.h>
#include <pacbio/exception/InvalidEvaluatorException.h>

using namespace std;

namespace PacBio {
namespace Consensus {

PolishConfig::PolishConfig(const size_t iterations, const size_t separation,
                           const size_t neighborhood)
    : MaximumIterations(iterations)
    , MutationSeparation(separation)
    , MutationNeighborhood(neighborhood)
{
}

RepeatConfig::RepeatConfig(const size_t repeatSize, const size_t elementCount,
                           const size_t iterations)
    : MaximumRepeatSize{repeatSize}
    , MinimumElementCount{elementCount}
    , MaximumIterations{iterations}
{
}

void Mutations(vector<Mutation>* muts, const AbstractIntegrator& ai, const size_t start,
               const size_t end)
{
    constexpr auto bases = "ACGT";
    constexpr size_t nbases = 4;

    if (start == end) return;

    char last = (start > 0) ? ai[start - 1] : '\0';

    for (size_t i = start; i < end; ++i) {
        const char curr = ai[i];

        // insertions come before deletion/substitutions at site i, their End()
        // is i < i + 1
        for (size_t j = 0; j < nbases; ++j) {
            // if it's not a homopolymer insertion, or it's the first base of
            // one..
            if (bases[j] != last) muts->emplace_back(Mutation::Insertion(i, bases[j]));
        }

        // if we're the first in the homopolymer, we can delete
        if (curr != last) muts->emplace_back(Mutation::Deletion(i, 1));

        for (size_t j = 0; j < nbases; ++j) {
            if (bases[j] != curr) muts->emplace_back(Mutation::Substitution(i, bases[j]));
        }

        last = curr;
    }

    // if we are at the end, make sure we're not performing a terminal
    // homopolymer insertion
    for (size_t j = 0; j < nbases; ++j)
        if (bases[j] != last) muts->emplace_back(Mutation::Insertion(end, bases[j]));
}

vector<Mutation> Mutations(const AbstractIntegrator& ai, const size_t start, const size_t end)
{
    vector<Mutation> muts;
    Mutations(&muts, ai, start, end);
    return muts;
}

vector<Mutation> Mutations(const AbstractIntegrator& ai)
{
    return Mutations(ai, 0, ai.TemplateLength());
}

void RepeatMutations(vector<Mutation>* muts, const AbstractIntegrator& ai, const RepeatConfig& cfg,
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

vector<Mutation> RepeatMutations(const AbstractIntegrator& ai, const RepeatConfig& cfg,
                                 const size_t start, const size_t end)
{
    vector<Mutation> muts;
    RepeatMutations(&muts, ai, cfg, start, end);
    return muts;
}

vector<Mutation> RepeatMutations(const AbstractIntegrator& ai, const RepeatConfig& cfg)
{
    return RepeatMutations(ai, cfg, 0, ai.TemplateLength());
}

vector<Mutation> BestMutations(list<ScoredMutation>* scoredMuts, const size_t separation)
{
    vector<Mutation> result;

    // TODO handle 0-separation correctly
    if (separation == 0) throw invalid_argument("nonzero separation required");

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
                                 const AbstractIntegrator& ai, const size_t neighborhood)
{
    const size_t len = ai.TemplateLength();
    const auto clamp = [len](const int i) { return max(0, min<int>(len, i)); };

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
        Mutations(&result, ai, range.first, range.second);

    return result;
}

PolishResult Polish(AbstractIntegrator* const ai, const PolishConfig& cfg)
{
    vector<Mutation> muts = Mutations(*ai);
    hash<string> hashFn;
    size_t oldTpl = hashFn(*ai);
    set<size_t> history = {oldTpl};

    PolishResult result;

    for (size_t i = 0; i < cfg.MaximumIterations; ++i) {
        // find the best mutations given our parameters
        {
            list<ScoredMutation> scoredMuts;
            int mutationsTested = 0;
            bool hasNewInvalidEvaluator;

            // Compute new sets of possible mutations until no Evaluators are
            // being invalided.
            do {
                // Compute the LL only with the active Evaluators
                const double LL = ai->LL();

                hasNewInvalidEvaluator = false;
                try {
                    // Get set of possible mutations
                    for (const auto& mut : muts) {
                        const double ll = ai->LL(mut);
                        if (ll > LL) scoredMuts.emplace_back(mut.WithScore(ll));
                        ++mutationsTested;
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
            return result;
        }

        const size_t newTpl = hashFn(ApplyMutations(*ai, &muts));

        const auto diagnostics = [&result](AbstractIntegrator* ai) {
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
            muts = NearbyMutations(&applied, &muts, *ai, cfg.MutationNeighborhood);
        } else {
            ai->ApplyMutations(&muts);
            oldTpl = newTpl;
            result.mutationsApplied += muts.size();

            diagnostics(ai);

            // get the mutations for the next round
            muts = NearbyMutations(&muts, &muts, *ai, cfg.MutationNeighborhood);
        }

        // keep track of which templates we've seen
        history.insert(oldTpl);
    }

    return result;
}

PolishResult PolishRepeats(AbstractIntegrator* const ai, const RepeatConfig& cfg)
{
    PolishResult result;

    const auto diagnostics = [&result](AbstractIntegrator* ai) {
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
        throw invalid_argument("invalid value: probability not in [0,1]");
    else if (probability == 0.0)
        probability = numeric_limits<double>::min();

    return static_cast<int>(round(-10.0 * log10(probability)));
}

inline int ScoreSumToQV(const double scoreSum)
{
    return ProbabilityToQV(1.0 - 1.0 / (1.0 + scoreSum));
}

}  // anonymous namespace

vector<int> ConsensusQualities(AbstractIntegrator& ai)
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

QualityValues ConsensusQVs(AbstractIntegrator& ai)
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
    return QualityValues{move(quals), move(delQVs), move(insQVs), move(subQVs)};
}

}  // namespace Consensus
}  // namespace PacBio
