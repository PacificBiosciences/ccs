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

#include <algorithm>
#include <cfloat>
#include <map>
#include <string>
#include <vector>
#include <boost/format.hpp>

#include <ConsensusCore/Checksum.hpp>
#include <ConsensusCore/Mutation.hpp>
#include <ConsensusCore/Sequence.hpp>
#include <ConsensusCore/Utils.hpp>
#include <ConsensusCore/Arrow/MutationScorer.hpp>
#include <ConsensusCore/Arrow/MultiReadMutationScorer.hpp>
#include <ConsensusCore/Arrow/ContextParameters.hpp>
#include <ConsensusCore/Arrow/TemplateParameterPair.hpp>
#include <ConsensusCore/Arrow/SimpleRecursor.hpp>

#define MIN_FAVORABLE_SCOREDIFF 0.04  // Chosen such that 0.49 = 1 / (1 + exp(minScoreDiff))

namespace ConsensusCore {
namespace Arrow {
    /**
     Could the mutation change the contents of the portion of the
      template that is mapped to the read?
     @param read A read mapped to some portion of the template
     @param mut A mutation affecting the template
     @returns true if the read scores the mutation.
     @exception <#throws#>
     */
    //
    //
    template<typename R>
    bool MultiReadMutationScorer<R>::ReadScoresMutation(const MappedArrowRead& read, const Mutation& mut) const
    {
        int ts = read.TemplateStart;
        int te = read.TemplateEnd;
        int ms = mut.Start();
        int me = mut.End();
        if (mut.IsInsertion())
            return (ts <= me && ms <= te);
        return (ts < me && ms < te);    // Intervals intersect
    }


    //
    /**
     // Logic for turning a mutation to the global template space to
     // one in the coordinates understood by each individual mutation
     // scorer.  This involves translation, complementation, and also
     // possible clipping, if the mutation is not wholly within the
     // mapped read.
     @param <#parameter#>
     @returns <#retval#>
     */
    template<typename R>
    Mutation
    MultiReadMutationScorer<R>::OrientedMutation(const MappedArrowRead& mr,
                              const Mutation& mut)
    {
        using std::min;
        using std::max;

        // Clip mutation to bounds of mapped read, so that overhanging
        // multibase changes are handled correctly
        Mutation cmut(INSERTION, 0, 0, "N");
        if (mut.End() - mut.Start() > 1)
        {
            int cs, ce;
            cs = max(mut.Start(), mr.TemplateStart);
            ce = min(mut.End(), mr.TemplateEnd);
            if (mut.IsSubstitution())
            {
                std::string cNewBases = mut.NewBases().substr(cs-mut.Start(), ce-cs);
                cmut = Mutation(mut.Type(), cs, ce, cNewBases);
            }
            else
            {
                cmut = Mutation(mut.Type(), cs, ce, mut.NewBases());
            }
        }
        else
        {
            cmut = mut;
        }

        // Now orient
        if (mr.Strand == FORWARD_STRAND)
        {
            return Mutation(cmut.Type(),
                            cmut.Start() - mr.TemplateStart,
                            cmut.End() - mr.TemplateStart,
                            cmut.NewBases());
        }
        else
        {
            // This is tricky business
            int end   = mr.TemplateEnd - cmut.Start();
            int start = mr.TemplateEnd - cmut.End();
            return Mutation(cmut.Type(), start, end, ReverseComplement(cmut.NewBases()));
        }
    }



    template<typename R>
    MultiReadMutationScorer<R>::MultiReadMutationScorer(const ArrowConfig& config,
                                                        std::string tpl)
        : arrConfig_(config),
          fwdTemplate_(tpl, config.CtxParams),
          revTemplate_(ReverseComplement(tpl), config.CtxParams),
          reads_()
    {
        DEBUG_ONLY(CheckInvariants());
        ArrowConfigTable::const_iterator it;

    }

#if FALSE
    template<typename R>
    MultiReadMutationScorer<R>::MultiReadMutationScorer(const MultiReadMutationScorer<R>& other)
        : arrConfig_(other.arrConfig_),
          fwdTemplate_(other.fwdTemplate_),
          revTemplate_(other.revTemplate_),
          reads_()
    {
        // Make a deep copy of the readsAndScorers
        foreach (const ReadStateType& read, reads_)
        {
            reads_.push_back(ReadStateType(read));
        }

        DEBUG_ONLY(CheckInvariants());
    }
#endif

    template<typename R>
    MultiReadMutationScorer<R>::~MultiReadMutationScorer()
    {}

    template<typename R>
    int
    MultiReadMutationScorer<R>::TemplateLength() const
    {
        return (int)fwdTemplate_.tpl.length();
    }

    template<typename R>
    int
    MultiReadMutationScorer<R>::NumReads() const
    {
        return (int)reads_.size();
    }
    
    template<typename R>
    void
    MultiReadMutationScorer<R>::DumpAlphas()
    {
        for (int i = 0; i < NumReads(); ++i)
        {
            if (reads_[i].IsActive)
                DumpMatrix(*(reads_[i].Scorer->Alpha()),
                           "Alpha" + std::to_string(i + 1) + ".csv");
        }
    }
    
    template<typename R>
    const MappedArrowRead*
    MultiReadMutationScorer<R>::Read(int readIdx) const
    {
        return reads_[readIdx].IsActive ? reads_[readIdx].Read : NULL;
    }

    template<typename R>
    std::string
    MultiReadMutationScorer<R>::Template(StrandEnum strand)
    {
        return (strand == FORWARD_STRAND ? fwdTemplate_.tpl : revTemplate_.tpl);
    }

    template<typename R>
    WrappedTemplateParameterPair
    MultiReadMutationScorer<R>::Template(StrandEnum strand,
                                         int templateStart,
                                         int templateEnd)
    {
        int len = templateEnd - templateStart;
        if (strand == FORWARD_STRAND)
        {
            return fwdTemplate_.GetSubSection(templateStart, len);
        }
        else
        {
            return revTemplate_.GetSubSection(TemplateLength() - templateEnd, len);
        }
    }

    template<typename R>
    void
    MultiReadMutationScorer<R>::ApplyMutations(const std::vector<Mutation>& mutations)
    {
        DEBUG_ONLY(CheckInvariants());
        std::vector<int> mtp = TargetToQueryPositions(mutations, fwdTemplate_.tpl);
        fwdTemplate_.ApplyRealMutations(mutations, arrConfig_.CtxParams);
        revTemplate_.Reset(TemplateParameterPair(ReverseComplement(fwdTemplate_.tpl), arrConfig_.CtxParams));

        foreach (ReadStateType& rs, reads_)
        {
            try {
                int newTemplateStart = mtp[rs.Read->TemplateStart];
                int newTemplateEnd   = mtp[rs.Read->TemplateEnd];

                // reads (even inactive reads) will have their mapping coords updated
                rs.Read->TemplateStart = newTemplateStart;
                rs.Read->TemplateEnd   = newTemplateEnd;

                if (rs.IsActive)
                {
                    rs.Scorer->Template(Template(rs.Read->Strand,
                                                 newTemplateStart,
                                                 newTemplateEnd));
                }
            }
            catch (AlphaBetaMismatchException& e)
            {
                rs.IsActive = false;
            }
        }
        DEBUG_ONLY(CheckInvariants());
    }
    /**
     Add a read to the multi read mutation scorer
     @param mr A mapped read to addd
     @param threshold A criteria to ensure that not so much memory is used. If the read consumes more than threshold * (I+1)*(J+1) memory while filling the forward and backwards matrices, the read is not added to the scorer.
     @returns If the read was added to the scorer.
     @exception <#throws#>
     */
    template<typename R>
    AddReadResult MultiReadMutationScorer<R>::AddRead(const MappedArrowRead& mr, double threshold)
    {
        AddReadResult res = SUCCESS;
        DEBUG_ONLY(CheckInvariants());
        RecursorType recursor(arrConfig_.MdlParams, mr, Template(mr.Strand, mr.TemplateStart, mr.TemplateEnd), arrConfig_.Banding);

        ScorerType* scorer;
        try
        {
            scorer = new MutationScorer<R>(recursor);
        }
        catch (AlphaBetaMismatchException& e)
        {
            scorer = NULL;
            res = ALPHABETAMISMATCH;
        }

        if (scorer != NULL && threshold < 1.0f)
        {
            int I = recursor.read_.Length();
            int J = recursor.tpl_.Length();
            int maxSize = static_cast<int>(0.5 + threshold * (I + 1) * (J + 1));
            
            if (scorer->Alpha()->AllocatedEntries() >= maxSize ||
                scorer->Beta()->AllocatedEntries() >= maxSize)
            {
                res = AddReadResult::MEM_FAIL;
                delete scorer;
                scorer = NULL;
            }
        }

        bool isActive = scorer != NULL;
        reads_.push_back(ReadStateType(new MappedArrowRead(mr), scorer, isActive));
        DEBUG_ONLY(CheckInvariants());
        return res;
    }
    /**
     Add a read to the scorer using the default memory threshold.
     @param mr A mapped read to add
     @returns If the read was included
     */
    template<typename R>
    AddReadResult MultiReadMutationScorer<R>::AddRead(const MappedArrowRead& mr)
    {
        DEBUG_ONLY(CheckInvariants());
        return AddRead(mr, arrConfig_.AddThreshold);
    }

    template<typename R>
    double MultiReadMutationScorer<R>::Score(const Mutation& m, double fastScoreThreshold)
    {
        // Apply virtual mutations
        // First to virtually apply the mutation to the underlying type
        fwdTemplate_.ApplyVirtualMutation(m, arrConfig_.CtxParams);
        // Now create and apply the reverse complement mutation.
        int end   = (int)fwdTemplate_.tpl.length() - m.Start(); //Used to be int end   = mr.TemplateEnd - cmut.Start();
        int start = (int)fwdTemplate_.tpl.length() - m.End();
        auto rc_m = Mutation(m.Type(), start, end, ReverseComplement(m.NewBases()));
        revTemplate_.ApplyVirtualMutation(rc_m, arrConfig_.CtxParams);
    
        // Now score the mutation on each of them.
        double sum = 0;
        foreach (const ReadStateType& rs, reads_)
        {
            if (rs.IsActive && ReadScoresMutation(*rs.Read, m))
            {
                Mutation orientedMut = OrientedMutation(*rs.Read, m);
                sum += (rs.Scorer->ScoreMutation(orientedMut) -
                        rs.Scorer->Score());
            }
            if (sum < fastScoreThreshold)
                break;
        }
        // Clear the virtual mutation
        fwdTemplate_.ClearVirtualMutation();
        revTemplate_.ClearVirtualMutation();
        assert(!fwdTemplate_.VirtualMutationActive() && !revTemplate_.VirtualMutationActive());
        return sum;
    }

    template<typename R>
    double MultiReadMutationScorer<R>::Score(MutationType mutationType,
                                             int position, char base)
    {
        Mutation m(mutationType, position, base);
        return Score(m);
    }

    template<typename R>
    double MultiReadMutationScorer<R>::FastScore(const Mutation& m)
    {
        return Score(m, arrConfig_.FastScoreThreshold);
    }

    template<typename R>
    std::vector<double>
    MultiReadMutationScorer<R>::Scores(const Mutation& m, double unscoredValue)
    {
        std::vector<double> scoreByRead;
        foreach (const ReadStateType& rs, reads_)
        {
            if (rs.IsActive && ReadScoresMutation(*rs.Read, m))
            {
                Mutation orientedMut = OrientedMutation(*rs.Read, m);
                scoreByRead.push_back(rs.Scorer->ScoreMutation(orientedMut) -
                                      rs.Scorer->Score());
            }
            else
            {
                scoreByRead.push_back(unscoredValue);
            }
        }
        return scoreByRead;
    }

    template<typename R>
    std::vector<double> MultiReadMutationScorer<R>::Scores(MutationType mutationType,
                                                          int position, char base,
                                                          double unscoredValue)
    {
        Mutation m(mutationType, position, base);
        return Scores(m, unscoredValue);
    }

    template<typename R>
    bool MultiReadMutationScorer<R>::IsFavorable(const Mutation& m)
    {
        double sum = Score(m);
        return (sum > MIN_FAVORABLE_SCOREDIFF);
    }

    template<typename R>
    bool MultiReadMutationScorer<R>::FastIsFavorable(const Mutation& m)
    {
        double sum = FastScore(m);
        return (sum > MIN_FAVORABLE_SCOREDIFF);
    }


    template<typename R>
    std::vector<int> MultiReadMutationScorer<R>::AllocatedMatrixEntries() const
    {
        std::vector<int> allocatedCounts;
        for (int i = 0; i < (int)reads_.size(); i++)
        {
            int n = AlphaMatrix(i)->AllocatedEntries() + BetaMatrix(i)->AllocatedEntries() ;
            allocatedCounts.push_back(n);
        }
        return allocatedCounts;
    }


    template<typename R>
    std::vector<int> MultiReadMutationScorer<R>::UsedMatrixEntries() const
    {
        std::vector<int> usedCounts;
        for (int i = 0; i < (int)reads_.size(); i++)
        {
            int n =  AlphaMatrix(i)->UsedEntries() + BetaMatrix(i)->UsedEntries();
            usedCounts.push_back(n);
        }
        return usedCounts;
    }


    template<typename R>
    const AbstractMatrix<double>* MultiReadMutationScorer<R>::AlphaMatrix(int i) const
    {
        return reads_[i].Scorer->Alpha();
    }


    template<typename R>
    const AbstractMatrix<double>* MultiReadMutationScorer<R>::BetaMatrix(int i) const
    {
        return reads_[i].Scorer->Beta();
    }


    template<typename R>
    std::vector<int> MultiReadMutationScorer<R>::NumFlipFlops() const
    {
        std::vector<int> nFlipFlops;
        foreach (const ReadStateType& rs, reads_)
        {
            nFlipFlops.push_back(rs.Scorer->NumFlipFlops());
        }
        return nFlipFlops;
    }


    template<typename R>
    double MultiReadMutationScorer<R>::BaselineScore() const
    {
        double sum = 0;
        foreach (const ReadStateType& rs, reads_)
        {
            if (rs.IsActive) sum += rs.Scorer->Score();
        }
        return sum;
    }


    template<typename R>
    std::vector<double> MultiReadMutationScorer<R>::BaselineScores() const
    {
        std::vector<double> scoreByRead;
        foreach (const ReadStateType& rs, reads_)
        {
            if (rs.IsActive) scoreByRead.push_back(rs.Scorer->Score());
        }
        return scoreByRead;
    }
    
    template<typename R>
    void MultiReadMutationScorer<R>::DumpMatrix(const ScaledSparseMatrixD& mat, const std::string& fname)
    {
        if (mat.Rows() == 0 || mat.Columns() == 0) return;
        std::ofstream myfile;
        myfile.open(fname);
        for(int i=0; i< mat.Rows(); i++)
        {
            myfile << mat(i,0);
            for(int j=1; j<mat.Columns(); j++)
            {
                myfile << "," << mat(i,j);
            }
            myfile << std::endl;
        }
        myfile << mat.GetLogScale(0);
        for (int j=1; j < mat.Columns(); j++)
        {
            myfile << "," << mat.GetLogScale(j);
        }
        myfile << std::endl;
        myfile.close();
    }


    template<typename R>
    void MultiReadMutationScorer<R>::CheckInvariants() const
    {
#ifndef NDEBUG
        assert(!fwdTemplate_.VirtualMutationActive() && !revTemplate_.VirtualMutationActive());
        assert(revTemplate_.tpl == ReverseComplement(fwdTemplate_.tpl));
        foreach (const ReadStateType& rs, reads_)
        {
            rs.CheckInvariants();
            if (rs.IsActive) {
                //assert(rs.Scorer->Template().tpl == Template(rs.Read->Strand,
                //                                         rs.Read->TemplateStart,
                //                                         rs.Read->TemplateEnd).tpl);
                assert(0 <= rs.Read->TemplateStart && rs.Read->TemplateStart <= fwdTemplate_.tpl.size());
                assert(0 <= rs.Read->TemplateEnd && rs.Read->TemplateEnd <= fwdTemplate_.tpl.size());
                assert(rs.Read->TemplateStart <= rs.Read->TemplateEnd);
            }
        }
#endif  // !NDEBUG
    }


    template<typename R>
    std::string MultiReadMutationScorer<R>::ToString() const
    {
        std::stringstream ss;

        ss << "Template: " << fwdTemplate_.tpl << std::endl;
        ss << "Score: " << BaselineScore() << std::endl;

        ss << "Reads:" << std::endl;
        foreach (const ReadStateType& rs, reads_)
        {
            ss << "\t" << rs.ToString() << std::endl;
        }
        return ss.str();
    }


    namespace detail {

        template<typename ScorerType>
        ReadState<ScorerType>::ReadState(MappedArrowRead* read,
                                         ScorerType* scorer,
                                         bool isActive)
            : Read(read),
              Scorer(scorer),
              IsActive(isActive)
        {
            CheckInvariants();
        }

        template<typename ScorerType>
        ReadState<ScorerType>::ReadState(const ReadState& other)
            : Read(NULL),
              Scorer(NULL),
              IsActive(other.IsActive)
        {
            if (other.Read != NULL) Read = new MappedArrowRead(*other.Read);
            if (other.Scorer != NULL) Scorer = new ScorerType(*other.Scorer);
            CheckInvariants();
        }

        template<typename ScorerType>
        ReadState<ScorerType>::~ReadState()
        {
            if (Read != NULL) delete Read;
            if (Scorer != NULL) delete Scorer;
        }

        template<typename ScorerType>
        void ReadState<ScorerType>::CheckInvariants() const
        {
#ifndef NDEBUG
            if (IsActive)
            {
                assert(Read != NULL && Scorer != NULL);
                assert((int)Scorer->Template().Length() ==
                       Read->TemplateEnd - Read->TemplateStart);
            }
#endif  // !NDEBUG
        }

        template<typename ScorerType>
        std::string ReadState<ScorerType>::ToString() const
        {
            std::string score;
            if (IsActive)
            {
                score = (boost::format(" (Score= %0.2f)") % Scorer->Score()).str();
            }
            else
            {
                score = "*INACTIVE*";
            }
            return Read->ToString() + score;

        }
    }
    
    template class MultiReadMutationScorer<ArrowRecursor>;
}
}
