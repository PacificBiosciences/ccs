// Copyright (c) 2011-2014, Pacific Biosciences of California, Inc.
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

// Author: David Alexander, Lance Hepler

#include <ConsensusCore/MutationEnumerator.hpp>
#include <ConsensusCore/Mutation.hpp>
#include <ConsensusCore/Utils.hpp>
#include <ConsensusCore/Logging.hpp>

#include <algorithm>
#include <boost/functional/hash.hpp>
#include <boost/tuple/tuple.hpp>
#include <cmath>
#include <limits>
#include <set>
#include <string>
#include <utility>
#include <vector>


namespace ConsensusCore
{
    namespace {  // PRIVATE
    struct RefineRepeatOptions : RefineOptions
    {
        explicit RefineRepeatOptions(int repeatLength, int minRepeatElements)
            : RepeatLength(repeatLength)
            , MinRepeatElements(minRepeatElements)
        {
            MaximumIterations = 1;
        }

        int RepeatLength;
        int MinRepeatElements;
    };

    std::vector<ScoredMutation>
    DeleteRange(std::vector<ScoredMutation> input, int rStart, int rEnd)
    {
        std::vector<ScoredMutation> output;
        foreach (ScoredMutation s, input)
        {
            int pos = s.Start();
            if (!(rStart <= pos && pos <= rEnd))
            {
                output.push_back(s);
            }
        }
        return output;
    }

    bool ScoreComparer(const ScoredMutation& i, const ScoredMutation& j)
    {
        return i.Score() < j.Score();
    }

    //    Given a list of (mutation, score) tuples, this utility method
    //    greedily chooses the highest scoring well-separated elements.  We
    //    use this to avoid applying adjacent high scoring mutations, which
    //    are the rule, not the exception.  We only apply the best scoring one
    //    in each neighborhood, and then revisit the neighborhoods after
    //    applying the mutations.
    //
    //    This is highly unoptimized.  It is not in the critical path.
    std::vector<ScoredMutation>
    BestSubset(std::vector<ScoredMutation> input, int mutationSeparation)
    {
        if (mutationSeparation == 0)
            return input;

        std::vector<ScoredMutation> output;

        while (!input.empty())
        {
            ScoredMutation& best = *std::max_element(input.begin(), input.end(), ScoreComparer);
            output.push_back(best);
            int nStart = best.Start() - mutationSeparation;
            int nEnd = best.Start() + mutationSeparation;
            input = DeleteRange(input, nStart, nEnd);
        }

        return output;
    }


    // Sadly and annoyingly there is no covariance on std::vector in C++, so we have
    // to explicitly project back down to the superclass type to use the APIs as written.
    std::vector<Mutation>
    ProjectDown(const std::vector<ScoredMutation>& smuts)
    {
        return std::vector<Mutation>(smuts.begin(), smuts.end());
    }


    int ProbabilityToQV(double probability)
    {
        if (probability < 0.0 || probability > 1.0)
            throw std::invalid_argument("invalid value: probability not in [0,1]");
        else if (probability == 0.0)
            probability = std::numeric_limits<double>::min();

        return static_cast<int>(round(-10.0 * log10(probability)));
    }

    template <typename E, typename O>
    E MutationEnumerator(const std::string& tpl, const O& opts)
    {
        return E(tpl);
    }

    //
    // this MUST go last to properly specialize the MutationEnumerator
    //
    template <>
    RepeatMutationEnumerator
    MutationEnumerator<>(const std::string& tpl, const RefineRepeatOptions& opts)
    {
        return RepeatMutationEnumerator(
            tpl,
            opts.RepeatLength,
            opts.MinRepeatElements);
    }

    template <typename E, typename M, typename O>
    bool AbstractRefineConsensus(M& mms, const O& opts)
    {
        bool isConverged = false;
        float score = mms.BaselineScore();
        boost::hash<std::string> hash;
        std::set<size_t> tplHistory;

        std::vector<ScoredMutation> favorableMutsAndScores;

        for (int iter = 0; iter < opts.MaximumIterations; iter++)
        {
            LDEBUG << "Round " << iter;
            LDEBUG << "State of MMS: " << std::endl << mms.ToString();

            if (tplHistory.find(hash(mms.Template())) != tplHistory.end())
            {
                LDEBUG << "Cycle detected!";
            }

            if (mms.BaselineScore() < score)
            {
                LDEBUG << "Score decrease";  // Usually recoverable
            }
            score = mms.BaselineScore();

            //
            // Try all mutations in iteration 0.  In subsequent iterations, try mutations
            // nearby those used in previous iteration.
            //
            E mutationEnumerator = MutationEnumerator<E, O>(mms.Template(), opts);
            std::vector<Mutation> mutationsToTry;
            if (iter == 0) {
                mutationsToTry = mutationEnumerator.Mutations();
            }
            else
            {
                mutationsToTry = UniqueNearbyMutations(mutationEnumerator,
                                                       ProjectDown(favorableMutsAndScores),
                                                       opts.MutationNeighborhood);
            }

            //
            // Screen for favorable mutations.  If none, we are done (converged).
            //
            favorableMutsAndScores.clear();
            foreach (const Mutation& m, mutationsToTry)
            {
                if (mms.FastIsFavorable(m)) {
                    float mutScore = mms.Score(m);
                    favorableMutsAndScores.push_back(m.WithScore(mutScore));
                }
            }
            if (favorableMutsAndScores.empty())
            {
                isConverged = true;
                break;
            }

            //
            // Go with the "best" subset of well-separated high scoring mutations
            //
            std::vector<ScoredMutation> bestSubset = BestSubset(favorableMutsAndScores,
                                                                opts.MutationSeparation);

            //
            // Attempt to avoid cycling.  We could do a better job here.
            //
            if (bestSubset.size() > 1)
            {
                std::string nextTpl = ApplyMutations(ProjectDown(bestSubset), mms.Template());
                if (tplHistory.find(hash(nextTpl)) != tplHistory.end())
                {
                    LDEBUG << "Attempting to avoid cycle";
                    bestSubset = std::vector<ScoredMutation>(bestSubset.begin(),
                                                             bestSubset.begin() + 1);
                }
            }

            LDEBUG << "Applying mutations:";
            foreach (const ScoredMutation& smut, bestSubset)
            {
                LDEBUG << "\t" << smut;
            }

            tplHistory.insert(hash(mms.Template()));
            mms.ApplyMutations(ProjectDown(bestSubset));
        }

        return isConverged;
    }
    }  // PRIVATE


    template<typename MultiReadScorerType>
    bool RefineConsensus(MultiReadScorerType& mms, const RefineOptions& opts)
    {
        return AbstractRefineConsensus<UniqueSingleBaseMutationEnumerator>(mms, opts);
    }


    template<typename MultiReadScorerType>
    void RefineRepeats
    (MultiReadScorerType& mms, int repeatLength, int minRepeatElements)
    {
        RefineRepeatOptions opts(repeatLength, minRepeatElements);
        AbstractRefineConsensus<RepeatMutationEnumerator>(mms, opts);
    }


    template<typename MultiReadScorerType>
    std::vector<int> ConsensusQVs(MultiReadScorerType& mms)
    {
        std::vector<int> QVs;
        UniqueSingleBaseMutationEnumerator mutationEnumerator(mms.Template());
        for (size_t pos = 0; pos < mms.Template().length(); pos++)
        {
            double scoreSum = 0.0;
            foreach (const Mutation& m, mutationEnumerator.Mutations(pos, pos + 1))
            {
                // TODO (lhepler): this is dumb, but untestable mutations,
                //   aka insertions at ends, cause all sorts of weird issues
                double score = mms.Score(m);
                if (score < 0.0)
                {
                    scoreSum += exp(score);
                }
            }
            QVs.push_back(ProbabilityToQV(1.0 - 1.0 / (1.0 + scoreSum)));
        }
        return QVs;
    }


#if 0
    Matrix<float> MutationScoresMatrix(mms)
    {
        NotYetImplemented();
    }


    Matrix<float> MutationScoresMatrix(mms, mutationsToScore)
    {
        NotYetImplemented();
    }
#endif
}
