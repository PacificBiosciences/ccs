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

// Author: David Alexander

#pragma once

#include <ConsensusCore/Types.hpp>
#include <ConsensusCore/Read.hpp>
#include <ConsensusCore/Matrix/AbstractMatrix.hpp>
#include <ConsensusCore/Quiver/MutationScorer.hpp>
#include <ConsensusCore/Quiver/QuiverConfig.hpp>
#include <ConsensusCore/Quiver/SseRecursor.hpp>

#include <boost/noncopyable.hpp>
#include <string>
#include <utility>
#include <vector>
#include <map>

namespace ConsensusCore {

    class AbstractMultiReadMutationScorer
    {
    protected:
        AbstractMultiReadMutationScorer() {}
        virtual ~AbstractMultiReadMutationScorer() {}

    public:
        virtual int TemplateLength() const = 0;
        virtual int NumReads() const = 0;
        virtual const MappedQvRead* Read(int readIndex) const = 0;

        virtual std::string Template(StrandEnum strand = FORWARD_STRAND) const = 0;
        virtual std::string Template(StrandEnum strand,
                                     int templateStart,
                                     int templateEnd) const = 0;

        virtual void ApplyMutations(const std::vector<Mutation>& mutations) = 0;

        // Reads provided must be clipped to the reference/scaffold window implied by the
        // template, however they need not span the window entirely---nonspanning reads
        // must be provided with (0-based) template start/end coordinates.
        virtual bool AddRead(const MappedQvRead& mappedRead, float threshold) = 0;
        virtual bool AddRead(const MappedQvRead& mappedRead) = 0;

        virtual float Score(const Mutation& m) const = 0;
        virtual float FastScore(const Mutation& m) const = 0;

        // Return a vector (of length NumReads) of the difference in
        // the score of each read caused by the template mutation.  In
        // the case where the mutation cannot be scored for a read
        // (i.e., it is too close to the end of the template, or the
        // read does not span the mutation site) that entry in the
        // vector is 0
        virtual std::vector<float> Scores(const Mutation& m, float unscoredValue) const = 0;
        virtual std::vector<float> Scores(const Mutation& m) const = 0;

        virtual bool IsFavorable(const Mutation& m) const = 0;
        virtual bool FastIsFavorable(const Mutation& m) const = 0;

        // Rough estimate of memory consumption of scoring machinery
        virtual std::vector<int> AllocatedMatrixEntries() const = 0;
        virtual std::vector<int> UsedMatrixEntries() const = 0;
        virtual const AbstractMatrixF* AlphaMatrix(int i) const = 0;
        virtual const AbstractMatrixF* BetaMatrix(int i) const = 0;
        virtual std::vector<int> NumFlipFlops() const = 0;

#if !defined(SWIG) || defined(SWIGCSHARP)
        // Alternate entry points for C# code, not requiring zillions of object
        // allocations.
        virtual float Score(MutationType mutationType,
                            int position,
                            const std::string& newBases) const = 0;
        virtual std::vector<float> Scores(MutationType mutationType,
                                          int position,
                                          const std::string& newBases,
                                          float unscoredValue) const = 0;
        virtual std::vector<float> Scores(MutationType mutationType,
                                          int position,
                                          const std::string& newBases) const = 0;
#endif

        // Return the actual sum of scores for the current template.
        // TODO(dalexander): need to refactor to make the semantics of
        // the various "Score" functions clearer.
        virtual float BaselineScore() const = 0;
        virtual std::vector<float> BaselineScores() const = 0;


        virtual std::string ToString() const = 0;
    };


    bool ReadScoresMutation(const MappedQvRead& mr, const Mutation& mut);
    Mutation OrientedMutation(const MappedQvRead& mr, const Mutation& mut);


    namespace detail {
        template<typename ScorerType>
        struct ReadState
        {
            MappedQvRead* Read;
            ScorerType* Scorer;
            bool IsActive;

            ReadState(MappedQvRead* read,
                      ScorerType* scorer,
                      bool isActive);

            ReadState(const ReadState& other);
            ~ReadState();
            void CheckInvariants() const;
            std::string ToString() const;
        };
    }

    template<typename R>
    class MultiReadMutationScorer : public AbstractMultiReadMutationScorer
    {
    public:
        typedef R                                         RecursorType;
        typedef typename R::EvaluatorType                 EvaluatorType;
        typedef typename ConsensusCore::MutationScorer<R> ScorerType;
        typedef typename detail::ReadState<ScorerType>    ReadStateType;

    public:
        MultiReadMutationScorer(const QuiverConfigTable& paramsByChemistry, std::string tpl);
        MultiReadMutationScorer(const MultiReadMutationScorer<R>& scorer);
        virtual ~MultiReadMutationScorer();

        int TemplateLength() const;
        int NumReads() const;
        const MappedQvRead* Read(int readIndex) const;

        std::string Template(StrandEnum strand = FORWARD_STRAND) const;
        std::string Template(StrandEnum strand, int templateStart, int templateEnd) const;
        void ApplyMutations(const std::vector<Mutation>& mutations);

        // Reads provided must be clipped to the reference/scaffold window implied by the
        // template, however they need not span the window entirely---nonspanning reads
        // must be provided with (0-based) template start/end coordinates.
        bool AddRead(const MappedQvRead& mappedRead, float threshold);
        bool AddRead(const MappedQvRead& mappedRead);

        float Score(const Mutation& m) const;
        float FastScore(const Mutation& m) const;

        // Return a vector (of length NumReads) of the difference in
        // the score of each read caused by the template mutation.  In
        // the case where the mutation cannot be scored for a read
        // (i.e., it is too close to the end of the template, or the
        // read does not span the mutation site) that entry in the
        // vector is -FLT_MAX, which is to be interpreted as NA.
        std::vector<float> Scores(const Mutation& m, float unscoredValue) const;
        std::vector<float> Scores(const Mutation& m) const
        {
            return Scores(m, 0.0f);
        }

        bool IsFavorable(const Mutation& m) const;
        bool FastIsFavorable(const Mutation& m) const;

        // Rough estimate of memory consumption of scoring machinery
        std::vector<int> AllocatedMatrixEntries() const;
        std::vector<int> UsedMatrixEntries() const;
        const AbstractMatrixF* AlphaMatrix(int i) const;
        const AbstractMatrixF* BetaMatrix(int i) const;
        std::vector<int> NumFlipFlops() const;

#if !defined(SWIG) || defined(SWIGCSHARP)
        // Alternate entry points for C# code, not requiring zillions of object
        // allocations.
        float Score(MutationType mutationType,
                    int position,
                    const std::string& newBases) const;
        std::vector<float> Scores(MutationType mutationType,
                                  int position,
                                  const std::string& newBases,
                                  float unscoredValue) const;
        std::vector<float> Scores(MutationType mutationType,
                                  int position,
                                  const std::string& newBases) const
        {
            return Scores(mutationType, position, newBases, 0.0f);
        }
#endif

    public:
        // Return the actual sum of scores for the current template.
        // TODO(dalexander): need to refactor to make the semantics of
        // the various "Score" functions clearer.
        float BaselineScore() const;
        std::vector<float> BaselineScores() const;

    public:
        std::string ToString() const;

    private:
        void CheckInvariants() const;

    private:
        QuiverConfigTable quiverConfigByChemistry_;
        float fastScoreThreshold_;
        std::string fwdTemplate_;
        std::string revTemplate_;
        std::vector<ReadStateType> reads_;
    };

    typedef MultiReadMutationScorer<SparseSseQvRecursor> \
      SparseSseQvMultiReadMutationScorer;
    typedef MultiReadMutationScorer<SparseSseQvSumProductRecursor> \
      SparseSseQvSumProductMultiReadMutationScorer;
}
