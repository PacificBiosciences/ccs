// Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.
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

#pragma once

#include <boost/noncopyable.hpp>
#include <string>
#include <utility>
#include <vector>
#include <map>
#include <cfloat>

#include <ConsensusCore/Read.hpp>
#include <ConsensusCore/Types.hpp>
#include <ConsensusCore/Matrix/AbstractMatrix.hpp>
#include <ConsensusCore/Matrix/ScaledMatrix.hpp>
#include <ConsensusCore/Matrix/SparseMatrix.hpp>
#include <ConsensusCore/Arrow/MutationScorer.hpp>
#include <ConsensusCore/Arrow/ArrowConfig.hpp>
#include <ConsensusCore/Arrow/detail/Combiner.hpp>

namespace ConsensusCore {
namespace Arrow {

    enum AddReadResult {SUCCESS, ALPHABETAMISMATCH, MEM_FAIL, OTHER};
    static const char* AddReadResultNames[] = { "SUCCESS", "ALPHA/BETA MISMATCH", "EXCESSIVE MEMORY USAGE", "OTHER" };

    namespace detail {
        template<typename ScorerType>
        struct ReadState
        {
            MappedArrowRead* Read;
            ScorerType* Scorer;
            bool IsActive;

            ReadState(MappedArrowRead* read,
                      ScorerType* scorer,
                      bool isActive);

            ReadState(const ReadState& other);
            ~ReadState();
            void CheckInvariants() const;
            std::string ToString() const;
        };
    }

    template<typename R>
    class MultiReadMutationScorer
    {

        public:
            typedef R                                                RecursorType;
            typedef typename ConsensusCore::Arrow::MutationScorer<R> ScorerType;
            typedef typename detail::ReadState<ScorerType>           ReadStateType;

        public:
            MultiReadMutationScorer(const ArrowConfig& config, std::string tpl);
            // Move constructor and assignment operator
            MultiReadMutationScorer(MultiReadMutationScorer&& rhs) = default;
            MultiReadMutationScorer& operator=(MultiReadMutationScorer&& rhs) = default;
#if FALSE
            MultiReadMutationScorer(const MultiReadMutationScorer<R>& scorer);
#endif
            virtual ~MultiReadMutationScorer();

            bool ReadScoresMutation(const MappedArrowRead& mr, const Mutation& mut) const;
            Mutation OrientedMutation(const MappedArrowRead& mr, const Mutation& mut);

            int TemplateLength() const;
            int NumReads() const;
            const MappedArrowRead* Read(int readIndex) const;

            std::string Template(StrandEnum strand = FORWARD_STRAND);

            void ApplyMutations(const std::vector<Mutation>& mutations);

            // Reads provided must be clipped to the reference/scaffold window implied by the
            // template, however they need not span the window entirely---nonspanning reads
            // must be provided with (0-based) template start/end coordinates.
            AddReadResult AddRead(const MappedArrowRead& mappedRead, double threshold);
            AddReadResult AddRead(const MappedArrowRead& mappedRead);

            double Score(const Mutation& m, double scoreThreshold = -DBL_MAX);
            double FastScore(const Mutation& m);

            // Return a vector (of length NumReads) of the difference in
            // the score of each read caused by the template mutation.  In
            // the case where the mutation cannot be scored for a read
            // (i.e., it is too close to the end of the template, or the
            // read does not span the mutation site) that entry in the
            // vector is -FLT_MAX, which is to be interpreted as NA.
            std::vector<double> Scores(const Mutation& m, double unscoredValue);
            std::vector<double> Scores(const Mutation& m)
            {
                return Scores(m, 0.0);
            }



            bool IsFavorable(const Mutation& m);
            bool FastIsFavorable(const Mutation& m);

            // Rough estimate of memory consumption of scoring machinery
            std::vector<int> AllocatedMatrixEntries() const;
            std::vector<int> UsedMatrixEntries() const;
            const AbstractMatrixD* AlphaMatrix(int i) const;
            const AbstractMatrixD* BetaMatrix(int i) const;
            std::vector<int> NumFlipFlops() const;

    #if !defined(SWIG) || defined(SWIGCSHARP)
            // Alternate entry points for C# code, not requiring zillions of object
            // allocations.
            double Score(MutationType mutationType, int position, char base);

            std::vector<double> Scores(MutationType mutationType,
                                      int position, char base,
                                      double unscoredValue);

            std::vector<double> Scores(MutationType mutationType,
                                      int position, char base)
            {
                return Scores(mutationType, position, base, 0.0);
            }
    #endif

        public:
            // Return the actual sum of scores for the current template.
            // TODO(dalexander): need to refactor to make the semantics of
            // the various "Score" functions clearer.
            double BaselineScore() const;
            std::vector<double> BaselineScores() const;
            void DumpAlphas();

        public:
            std::string ToString() const;

        private:
            void DumpMatrix(const ScaledSparseMatrixD& mat, const std::string& fname);
            void CheckInvariants() const;
            /*
             Create a thin wrapper around the base template with coordinates specific for this guy
             */
            WrappedTemplateParameterPair Template(StrandEnum strand, int templateStart, int templateEnd);

        protected:
            ArrowConfig arrConfig_;
            TemplateParameterPair fwdTemplate_;
            TemplateParameterPair revTemplate_;
            std::vector<ReadStateType> reads_;
        };

        // This is the main one that will be used in the loop.
        typedef MultiReadMutationScorer<ArrowRecursor> ArrowMultiReadMutationScorer;
}
}
