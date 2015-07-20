//
//  TemplateParameterPair.cpp
//  ConsensusCoreXCode
//
//  Created by Nigel Delaney on 3/10/15.
//  Copyright (c) 2015 Pacific Biosciences. All rights reserved.
//

#include <stdio.h>
#include <ConsensusCore/Arrow/TemplateParameterPair.hpp>

namespace ConsensusCore {
namespace Arrow {

    TemplateParameterPair::TemplateParameterPair(const std::string& tpl_,
                                                 const std::vector<TransitionParameters>& trans_probs_)
        : mutantPosition( NO_MUTATION_SET_FLAG )
        , mutationOffset( NO_MUTATION_SET_FLAG )
        , tpl(tpl_)
        , trans_probs(trans_probs_)
    {
        assert(tpl.size() == trans_probs.size());
    }

    TemplateParameterPair::TemplateParameterPair()
        : mutantPosition( NO_MUTATION_SET_FLAG )
        , mutationOffset( NO_MUTATION_SET_FLAG )
        , tpl()
        , trans_probs()
    {}

#if 0
    TemplateParameterPair::TemplateParameterPair(const TemplateParameterPair& other)
        : mutantPosition( NO_MUTATION_SET_FLAG )
        , mutationOffset( NO_MUTATION_SET_FLAG )
        , tpl(other.tpl)
        , trans_probs(other.trans_probs)
    {
        assert(tpl.size() == trans_probs.size());
    }
#endif

    TemplateParameterPair::TemplateParameterPair(const std::string& tpl_, const ContextParameters& ctx)
        : mutantPosition( NO_MUTATION_SET_FLAG )
        , mutationOffset( NO_MUTATION_SET_FLAG )
        , tpl(tpl_)
        , trans_probs(tpl_.size())
    {
        // Initialize the dinucleotide context values.
        for(int i = 0; i < (tpl.size() -1); i ++)
        {
            auto ps = ctx.GetParametersForContext(tpl.at(i), tpl.at(i+1));
            trans_probs[i] = ps;
        }
        // Filling in the end just to make the math easier.
        trans_probs[tpl.size()-1] = TransitionParameters();

        assert(tpl.size() == trans_probs.size());
    }

    void TemplateParameterPair::ClearVirtualMutation() {
        mutantPosition = NO_MUTATION_SET_FLAG;
        mutationOffset = NO_MUTATION_SET_FLAG;
        mutantBP[0] = '0';
        mutantBP[1] = '0';
        mutantParameters[0] = TransitionParameters();
        mutantParameters[1] = TransitionParameters();
    }

    void TemplateParameterPair::ApplyVirtualMutation(const Mutation& mut, const ContextParameters& ctx_params)
    {
        // When applying a mutation, we need to:
        // 1 - Update the arrays that hold the "temporary" contexts
        // 2 - Update the parameters used for indexing.

        // Clear old memory to aid bug detection.
        ClearVirtualMutation();

        int start = mut.Start();
        mutantPosition = start;
        assert(start >= 0 && start < (tpl.length() + 1));

        if (mut.IsSubstitution()) {
            assert(mut.NewBases().length() == 1);
            mutationOffset = 0;
            char newBP = mut.NewBases()[0];
            mutantBP[1] = newBP;
            if (start > 0)
            {
                mutantBP[0] = this->tpl[start - 1];
                mutantParameters[0] = ctx_params.GetParametersForContext(tpl.at(start - 1), newBP);
            }
            if ((start + 1) < tpl.length()) {
                mutantParameters[1] = ctx_params.GetParametersForContext(newBP, tpl.at(start+1));
            }
        }
        else if (mut.IsDeletion()) {
            assert( (mut.End() - mut.Start()) == 1);
            mutationOffset = 1;
            auto org_length = tpl.length() - 1;
            // Three cases, at start, at end, and in middle.
            // If in middle, we update the prior position and ignore the removed position
            // If at the start, we only remove that position
            // If at the end, we remove the prior position

            if (start > 0 && start < org_length) {
                char nextBP = tpl.at(start + 1);
                char prevBP = tpl.at(start - 1);
                mutantBP[0] = prevBP;
                mutantBP[1] = nextBP;
                mutantParameters[0] = ctx_params.GetParametersForContext(prevBP, nextBP);
                mutantParameters[1] = trans_probs[start + 1];
            } else if (start == 0) { // At the start
                 char nextBP = tpl.at(start + 1);
                 mutantBP[1] = nextBP;
                 mutantParameters[1] = trans_probs[start + 1];
            } else if (start == org_length) { // At the end, parameters are erased.
                char prevBP = tpl.at(start - 1);
                mutantBP[0] = prevBP;
            }
        }
        else if (mut.IsInsertion()) {
            // Insertions indicate the position to place the base in (so if at X, what was at X is now at X+1).
            assert(mut.NewBases().length() == 1);
            mutationOffset = - 1;
            char newBP = mut.NewBases()[0];
            mutantBP[1] = newBP;
            // Need to update two parameters, the ones for this base and the one
            // before this base.  If inserted at the start, there is no base before
            if (start > 0) {
                char prevBP = tpl.at(start-1);
                mutantBP[0] = prevBP;
                mutantParameters[0] = ctx_params.GetParametersForContext(prevBP, newBP);
            }
            if (start < tpl.length()) {
                char oldBP = tpl.at(start);
                mutantParameters[1] = ctx_params.GetParametersForContext(newBP, oldBP);
            }
        }
    }




    // NOTE: Start is not just equal to Mut.Start() here because the location of a mutation can change
    // as earlier mutations are applied.

    // TODO: I am really not sure if this "apply in place" is any better than a simple reallocation and
    // refill, consider trying to do this a bit cleaner in one passes (how I do the rev comp).  The is legacy code.
    void
    TemplateParameterPair::_ApplyMutationInPlace(const Mutation& mut, int start, const ContextParameters& ctx_params)
    {
        if (mut.IsSubstitution())
        {
            tpl.replace(start, mut.End() - mut.Start(), mut.NewBases());
            if ((start + 1) < tpl.length()) {
                trans_probs[start] = ctx_params.GetParametersForContext(tpl.at(start), tpl.at(start+1));
            }
            if (start > 0) {
                trans_probs[start-1] = ctx_params.GetParametersForContext(tpl.at(start -1), tpl.at(start));
            }
        }
        else if (mut.IsDeletion())
        {
            assert( (mut.End() - mut.Start()) == 1);
            auto org_length = tpl.length() - 1;
            assert(start >=0 && start <= org_length);
            tpl.erase(start, mut.End() - mut.Start());
            // Three cases, at start, at end, and in middle.
            // If in middle, we update the prior position and delete the removed position
            // If at the start, we only remove that position
            // If at the end, we remove the prior position
            if (start > 0 && start < org_length) {
                trans_probs[start-1] = ctx_params.GetParametersForContext(tpl.at(start-1), tpl.at(start));
                trans_probs.erase(trans_probs.begin() + start, trans_probs.begin() + start + ( mut.End()- mut.Start()));
            } else if (start == 0) { // At the start
                trans_probs.erase(trans_probs.begin() + start, trans_probs.begin() + start + ( mut.End()- mut.Start()));
            } else if (start == org_length ) { // At the end
                trans_probs.erase(trans_probs.begin() + start - 1, trans_probs.begin() + start - 1  + ( mut.End()- mut.Start()));
            }
        }
        else if (mut.IsInsertion())
        {
            assert(start >= 0);
            // Add template base
            assert(tpl.size() == trans_probs.size());
            tpl.insert(start, mut.NewBases());
            if (start > trans_probs.size())
            {
                trans_probs.push_back(TransitionParameters());
            }
            else
            {
                trans_probs.insert(trans_probs.begin() + start, TransitionParameters());
            }
            assert(tpl.size() == trans_probs.size());
            // Need to update two parameters, the ones for this base and the one
            // before this base.  If inserted at the start, there is no base before
            if (start > 0) {
                trans_probs[start - 1] = ctx_params.GetParametersForContext(tpl.at(start-1), tpl.at(start));
            }
            // If inserted at the end, there is no "current" probabilities to update
            if (start < trans_probs.size()) {
                auto new_params = ctx_params.GetParametersForContext(tpl.at(start), tpl.at(start+1));
                trans_probs[start] = new_params;
            }
        }
    }

    void
    TemplateParameterPair::ApplyRealMutations(const std::vector<Mutation>& muts, const ContextParameters& ctx_params)
    {
        // Apply mutations
        std::vector<Mutation> sortedMuts(muts);
        std::sort(sortedMuts.begin(), sortedMuts.end());
        int runningLengthDiff = 0;
        foreach (const Mutation& mut, sortedMuts)
        {
            _ApplyMutationInPlace(mut, mut.Start() + runningLengthDiff, ctx_params);
            runningLengthDiff += mut.LengthDiff();
        }
    }

    WrappedTemplateParameterPair
    TemplateParameterPair::GetSubSection(int start, int len) {
        // Remove one as the transition parameters need to be one smaller.
        return WrappedTemplateParameterPair(this, start, len);
    }


    TemplateParameterPair
    TemplateParameterPair::GetReverseComplement(const ContextParameters&ctx) {
        //TODO: Verify return value optimization occurs.
        return TemplateParameterPair(ReverseComplement(tpl), ctx);

    }


    WrappedTemplateParameterPair::WrappedTemplateParameterPair(TemplateParameterPair* base, int start, int length) {
        this->base = base;
        this->start = start;
        this->length = length;
    }
}
}
