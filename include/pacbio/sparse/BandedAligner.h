// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
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

// Programmer: Brett Bowman

#pragma once

#include <algorithm>
#include <exception>

#include <seqan/seeds.h>
#include <seqan/journaled_set/score_biaffine.h>

namespace PacBio {
namespace SparseAlignment {

class BandedAligner
{
   
public:  // Derived types and values
   typedef seqan::Dna TAlphabet;
   typedef seqan::String<TAlphabet> TString;
   typedef seqan::StringSet<TString> TStringSet;

   typedef seqan::Seed<seqan::Simple> TSeed;
   typedef seqan::String<TSeed> TSeedChain;
   typedef std::pair<size_t, TSeedChain> THit;

   typedef seqan::AlignConfig<true, true, true, true> TAlignConfig;
   typedef seqan::Align<TString, seqan::ArrayGaps> TAlign;
   typedef std::pair<TAlign, seqan::AlignmentStats> TAlignPair;
   typedef seqan::Row<TAlign>::Type TRow;
   typedef seqan::Iterator<TRow>::Type TRowIter;
   typedef seqan::Score<long, seqan::BiAffine> TScoringScheme;

public:  // structors
    // Default constructor
    BandedAligner(const TStringSet& references,
                  const size_t flankingBuffer = 50,
                  const size_t minMatchLength = 5,
                  const size_t minQueryLength = 50,
                  const size_t bandExtension = 15)
        : references_(references)
        , flankingBuffer_(flankingBuffer)
        , minQueryLength_(minQueryLength)
        , minMatchLength_(minMatchLength)
        , bandExtension_(bandExtension)
    {
        using namespace seqan;
    
        // Sanity-check the two arguments with hard range requirements
        if (minMatchLength_ < 1)
            throw std::invalid_argument("minMatchLength must be >= 1");
        if (bandExtension_ < 1)
            throw std::invalid_argument("bandExtension must be >= 1");

        // NOTE: BLASR uses inverted scores, i.e. (-5, 6, 4, 5)
        //scoringScheme_ = TScoringScheme(5, -6, -4, -4, -5, -5);
    }

    // Move constructor
    BandedAligner(BandedAligner&& src) = default;
    // Copy constructor
    BandedAligner(const BandedAligner& src) = default;
    // Move assignment operator
    BandedAligner& operator=(BandedAligner&& rhs) = default;
    // Copy assignment operator
    BandedAligner& operator=(const BandedAligner& rhs) = default;
    // Destructor
    ~BandedAligner() = default;

public:  // non-modifying methods

    /// Using the seed-chain as a guide, align the banded region of probable 
    /// similarity between the query and the reference sequence and return 
    /// the alignment.
    ///
    /// Roughly equivalent to BLASR's SDPAlign
    ///
    /// \param  region  A seed for storing the final alignment region
    /// \param  query  The query sequence in which the seeds were found
    /// \param  hit  A pair with the reference sequence index and the
    ///              chain of seeds
    ///
    /// \return  TAlignPair  an Alignment/Stats pair object
    TAlignPair AlignHit(TSeed& region,
                        const TString& query,
                        const THit& hit)
    {
        using namespace seqan;

        // Initialize the object we will return
        TAlignPair pair;
        resize(rows(pair.first), 2);
    
        // Extract the input data for easier access
        const TString& reference = references_[hit.first];
        const TSeedChain& seedChain = hit.second;

        // Identify the window around the seeded region in the query to align
        size_t queryStart     = beginPositionH(front(seedChain));
        size_t queryMinStart  = ((queryStart <= flankingBuffer_) ? 0 : (queryStart - flankingBuffer_));
        size_t queryStartDiff = queryStart - queryMinStart;
        size_t queryEnd       = endPositionH(back(seedChain));
        size_t queryMaxEnd    = std::min(length(query), queryEnd + flankingBuffer_);
        size_t queryEndDiff   = queryMaxEnd - queryEnd;

        // Identify the window around the seeded region in the reference to align
        size_t refStart     = beginPositionV(front(seedChain));
        size_t refMinStart  = ((refStart <= flankingBuffer_) ? 0 : (refStart - flankingBuffer_));
        size_t refStartDiff = refStart - refMinStart;
        size_t refEnd       = endPositionV(back(seedChain));
        size_t refMaxEnd    = std::min(length(reference), refEnd + flankingBuffer_);
        size_t refEndDiff   = refMaxEnd - refEnd;

        // Given the size of the flanking regions around both sides of the anchors, 
        //    what is the maximum amount of buffer we can safely take from both?
        size_t startDiff = std::min(queryStartDiff, refStartDiff);
        size_t endDiff   = std::min(queryEndDiff, refEndDiff);

        // If the maximum size allow for the left-side flanks is different
        //    from what we previously calculated, recalculate
        size_t queryAlignStart = queryMinStart;
        if (startDiff != queryStartDiff)
            queryAlignStart += queryStartDiff - startDiff;

        size_t refAlignStart = refMinStart;
        if (startDiff != refStartDiff)
            refAlignStart += refStartDiff - startDiff;

        // Ditto the right-side flanks
        size_t queryAlignEnd = queryMaxEnd;
        if (endDiff != queryEndDiff)
            queryAlignEnd -= queryEndDiff - endDiff;

        size_t refAlignEnd = refMaxEnd;
        if (endDiff != refEndDiff)
            refAlignEnd -= refEndDiff - endDiff;

        // Extract infixes from sequences and set them as our alignment targets
        const auto queryInfix = infix(query, queryAlignStart, queryAlignEnd);
        const auto refInfix   = infix(reference, refAlignStart, refAlignEnd);
        assignSource(row(pair.first, 0), queryInfix);
        assignSource(row(pair.first, 1), refInfix);

        if (queryAlignStart == 0 && refAlignStart == 0)
        {
            // If both of our subsequences start at the beginning,
            //     we can use the seedChain as-is
            bandedChainAlignment(pair.first, seedChain, scoringScheme_, 
                    alignConfig_, bandExtension_);
        } else if (queryAlignStart > 10) {
            // Otherwise we need to left-shift the values in the seed chain
            const TSeedChain& shiftedChain = ShiftSeedChain(seedChain, 
                    queryAlignStart, refAlignStart);
            bandedChainAlignment(pair.first, shiftedChain, scoringScheme_, 
                    alignConfig_, bandExtension_);
        } else {
            // There appears to be a bug with bandedChainAlignment, where by it
            //    will fail to align with default values for the bandExtension_
            //    if the seed chain starts too close to the beginning of the
            //    query sequence.  This is a temporary hack to work around that.
            const TSeedChain& shiftedChain = ShiftSeedChain(seedChain, 
                    queryAlignStart, refAlignStart);
            bandedChainAlignment(pair.first, shiftedChain, scoringScheme_, 
                    alignConfig_, 10);
        }
       
        // Clip the alignment to the region with a good similarity
        ClipAlignment(&pair);
        
        // Now that the edges of the alignment are finalized we can calculate the region
        SetAlignmentRegion(pair.first, queryAlignStart, refAlignStart, &region);

        // Return the final alignment
        return pair;
    }

private:  // Private non-modifying functions

    /// The existing banded alignment algorithm provided by seqan crashes
    /// if given any overlapping seeds, while Lance's seed-chaining algorithm
    /// is both non-global and significantly faster.  In order to use the former
    /// with the later, it is first necessary to trim such overlaps from the data.
    ///
    /// Overlaps represent small repeat regions, usually homopolymers, whose
    /// edges are deliniated by the edges of the overlap.  So by swapping the
    /// overlapping start/end positions, we are effectively trimming the seeds
    /// back to the non-repetitive region and allowing the banded aligner to
    /// find the optimal path through the repeat.
    ///
    /// \param  chain  A seed-set of kmer matchs between the sequences
    void TrimOverlappingSeeds(TSeedChain* chain)
    {
        using namespace seqan;

        // SeqAn doesn't support "->at()", so we de-reference here
        TSeedChain& chainRef = *chain;

        // If each pair of adjacent seeds, I and J ...
        for (size_t j = 1; j < length(chainRef); ++j)
        {
            size_t i = j - 1;

            // ... if they over lap in either dimension ...
            if (endPositionH(chainRef[i]) > beginPositionH(chainRef[j]) ||
                endPositionV(chainRef[i]) > beginPositionV(chainRef[j]))
            {
                // ... swap their start and end positions in both dimensions
                size_t maxPosV = endPositionV(chainRef[i]);
                setEndPositionV(chainRef[i], beginPositionV(chainRef[j]));
                setBeginPositionV(chainRef[j], maxPosV);

                size_t maxPosH = endPositionH(chainRef[i]);
                setEndPositionH(chainRef[i], beginPositionH(chainRef[j]));
                setBeginPositionH(chainRef[j], maxPosH);
            }
        }
    }
    
    /// When identifying seeds from a seed chain that need to be removed
    /// for whatever reason, we don't want to do so immediately lest we
    /// confuse our indices.  Instead we set the horizonal (query) length
    /// to zero, to be removed afterward.
    ///
    /// \param  chain  A seed-set of kmer matchs between the sequences
    ///
    /// \return  int  The number of seeds so removed
    int RemoveZeroLengthSeeds(TSeedChain* chain)
    {
        using namespace seqan;

        // SeqAn doesn't support "->at()", so we de-reference here
        TSeedChain& chainRef = *chain;

        int size = length(chainRef);
        int curr = 0;

        for (int i = 0; i < size; ++i)
        {
            // If end is greater than beginning, its a still-valid seed.
            // Copy it to our current counter-position and increment 
            if (endPositionH(chainRef[i]) > beginPositionH(chainRef[i]))
            {
                chainRef[curr] = chainRef[i];
                curr++;
            }
        }

        // Erase the last "curr" seeds, which are now all duplicates
        for (int i = size - 1; i >= curr; --i)
            erase(chainRef, i);

        // Shrink the container to fit the new smaller dataset
        shrinkToFit(chainRef);
        return size - curr;
    }

    /// The existing banded alignment algorithm provided by seqan crashes
    /// if given any overlapping seeds, while Lance's seed-chaining algorithm
    /// is both non-global and significantly faster.  In order to use the former
    /// with the later, it is first necessary to remove any seeds entirely
    /// "contained" seeds from the seed chain.
    ///
    /// \param  chain  A seed-set of kmer matchs between the sequences
    ///
    /// \return  int  The number of seeds so removed
    int RemoveContainedSeeds(TSeedChain* chain)
    {
        using namespace seqan;
        
        // SeqAn doesn't support "->at()", so we de-reference here
        TSeedChain& chainRef = *chain;

        for (int m = length(chainRef) - 1; m > 0; --m)
        {
            int n = m - 1;
            bool mergeFound = false;
            while (mergeFound == false && n >= 0)
            {
                if ((beginPositionH(chainRef[n]) <  beginPositionH(chainRef[m]) &&
                    endPositionH(chainRef[n])   >= endPositionH(chainRef[m]))   ||
                    (beginPositionV(chainRef[n]) <  beginPositionV(chainRef[m]) &&
                    endPositionV(chainRef[n])   >= endPositionV(chainRef[m])))
                {
                    setEndPositionH(chainRef[m], beginPositionH(chainRef[m]));
                    mergeFound = true;
                }
                n--;
            }
        }

        // Remove zero-length seeds and return the count
        int numRemoved = RemoveZeroLengthSeeds(&chainRef);
        return numRemoved;
    }

    /// The existing banded chain alignment algorithm in SeqAn is global,
    /// whereas we want local alignments around our seed chains.  We can use 
    /// infixes of the sequences we want to work around this, but we have to
    /// adjust the start positions of our seeds accordingly so they are correct
    /// relative to the start of the infix.
    ///
    /// \param  input  A seed-set of kmer matchs between the sequences
    /// \param  leftShift  Number of horizontal bases to shift the seeds
    /// \param  upShift  Number of vertical bases to shift the seeds
    ///
    /// \return  TSeedChain  A copy of input with the starting locations
    TSeedChain ShiftSeedChain(const TSeedChain& input,
                              const size_t leftShift,
                              const size_t upShift)
    {
        using namespace seqan;
        TSeedChain output;

        for (size_t i = 0; i < length(input); ++i)
        {
            const TSeed& seed = input[i];
            const size_t queryStart = beginPositionH(seed) - leftShift;
            const size_t refStart   = beginPositionV(seed) - upShift;
            const size_t size       = seedSize(seed);
            appendValue(output, TSeed(queryStart, refStart, size));
        }

        return output;
    }

    /// Clip the alignment down to the highest accuracy core region
    /// longer than minQueryLength and update the stats object in
    /// the pair accordingly.
    ///
    /// \param  pair  An Alignment/Stats pair object
    ///
    /// \return  int  An error flag - 0 if successful, 1 if failure
    void ClipAlignment(TAlignPair* pair)
    {
        using namespace seqan;

        // Traverse the ends of the alignment for minimal match anchors
        TRow& queryRow = row(pair->first, 0);
        TRow& refRow   = row(pair->first, 1);
        int alignStart = FindLeftSideMinMatch(queryRow, refRow);                
        int alignEnd   = FindRightSideMinMatch(queryRow, refRow);   

        // If the alignment start or end positions don't make sense,
        //      this is probably a garbage alignment.  Compute the
        //      alignment stats as-is and abandon clipping.
        if (alignStart < 0 || alignEnd < 0)
        {
            computeAlignmentStats(pair->second, pair->first, scoringScheme_);
            return;
        }

        // Clip the alignment to the ends of the range and compute our stats
        setClippedBeginPosition(queryRow, alignStart);
        setClippedEndPosition(queryRow, alignEnd);
        setClippedBeginPosition(refRow, alignStart);
        setClippedEndPosition(refRow, alignEnd);
        computeAlignmentStats(pair->second, pair->first, scoringScheme_);

        // If the initial alignment is short, we succeeded but do not have room
        //     for refinement - exit here.
        if (LengthInSequence(queryRow, alignStart, alignEnd) <= minQueryLength_)
            return;

        // Otherwise we need to clear our clipping ahead of refinement
        clearClipping(queryRow);
        clearClipping(refRow);

        // Declare looping variables
        AlignmentStats newLeftStats, newRightStats;
        int newAlignStart, newAlignEnd;

        while (true) {
            // Clear the alignment stats from the previous iteration
            clear(newLeftStats);
            clear(newRightStats);

            // Calculate the next valid start and end position for the alignment
            newAlignStart = FindLeftSideMinMatch(queryRow, refRow, alignStart + minMatchLength_);                
            newAlignEnd   = FindRightSideMinMatch(queryRow, refRow, alignEnd - minMatchLength_);

            if (LengthInSequence(queryRow, newAlignStart, alignEnd) > minQueryLength_)
            {
                // If the alignment length in the query is long enough to be valid,
                //    calculate the accuracy of the new start with the old end
                setClippedBeginPosition(queryRow, newAlignStart);
                setClippedEndPosition(queryRow, alignEnd);
                setClippedBeginPosition(refRow, newAlignStart);
                setClippedEndPosition(refRow, alignEnd);
                computeAlignmentStats(newLeftStats, pair->first, scoringScheme_);
                clearClipping(queryRow);
                clearClipping(refRow);
            }

            if (LengthInSequence(queryRow, alignStart, newAlignEnd) > minQueryLength_)
            {
                // If the alignment length in the query is long enough to be valid,
                //    calculate the accuracy of the old start with the new end
                setClippedBeginPosition(queryRow, alignStart);
                setClippedEndPosition(queryRow, newAlignEnd);
                setClippedBeginPosition(refRow, alignStart);
                setClippedEndPosition(refRow, newAlignEnd);
                computeAlignmentStats(newRightStats, pair->first, scoringScheme_);
                clearClipping(queryRow);
                clearClipping(refRow);
            }

            if (newLeftStats.alignmentIdentity > newRightStats.alignmentIdentity &&
                newLeftStats.alignmentIdentity > pair->second.alignmentIdentity + minAccuracyImprovement_) {
                // If the new align-start produces a better alignment than the end,
                //     and the improvement is non-trivial, keep it
                alignStart = newAlignStart;
                pair->second = newLeftStats;
            } else if (newRightStats.alignmentIdentity > newLeftStats.alignmentIdentity &&
                       newRightStats.alignmentIdentity > pair->second.alignmentIdentity + minAccuracyImprovement_) {
                // If the new align-end produces a better alignment than the start,
                //     and the improvement is non-trivial, keep it
                alignEnd = newAlignEnd;
                pair->second = newRightStats;
            } else {
                // If both possible new positions are invalid or inferior, we 
                //    assume that the previous clipping was optimal
                setClippedBeginPosition(queryRow, alignStart);
                setClippedEndPosition(queryRow, newAlignEnd);
                setClippedBeginPosition(refRow, alignStart);
                setClippedEndPosition(refRow, newAlignEnd);
                break;
            }
        }

        return;
    }

    /// Compare two aligned sequences and, counting in from the left-end, find 
    /// the first position where they match for at least minMatchLength_ bases.
    ///
    /// \param  queryRow  The first row to compare
    /// \param  refRow  The second row to compare
    /// \param  startPos  Optional position to begin the search
    ///
    /// \return  int  The position in the row the match begins
    int FindLeftSideMinMatch(const TRow& queryRow,
                             const TRow& refRow,
                             const size_t startPos = 0)
    {
        size_t matchCount = 0;
        for (size_t i = startPos; i < seqan::length(queryRow); ++i)
        {
            // If either row is a gap, reset the counter 
            if (IsGap(refRow[i]) || IsGap(queryRow[i]))
            {
                matchCount = 0;
                continue;
            }

            // If the aligned characters aren't gaps but still don't match,
            //    it's a subsitution error
            if (queryRow[i] != refRow[i])
            {
                matchCount = 0;
                continue;
            }
            
            // If we made it this far it's a match - increment the counter
            matchCount++;

            // If we've seen enough matches in a row, return our position
            if (matchCount >= minMatchLength_)
            {
                return i - minMatchLength_ + 1;
            }
        }

        // If we made it this far, something went wrong - return an error
        return -1;
    }

    /// Compare two aligned sequences and, counting in from the right-end, find 
    /// the first position where they match for at least minMatchLength_ bases.
    ///
    /// \param  queryRow  The first row to compare
    /// \param  refRow  The second row to compare
    /// \param  startPos  Optional position to begin the search
    ///
    /// \return  int  The position in the alignment the match ends
    int FindRightSideMinMatch(const TRow& queryRow,
                              const TRow& refRow,
                              size_t startPos = 0)
    {
        // If not specified, default to starting at the right end
        if (startPos == 0)
            startPos = seqan::length(queryRow) - 1;

        size_t matchCount = 0;
        for (int i = startPos; i >= 0; --i)
        {
            // If either row is a gap, reset the counter 
            if (IsGap(refRow[i]) || IsGap(queryRow[i]))
            {
                matchCount = 0;
                continue;
            }

            // If the aligned characters aren't gaps but still don't match,
            //    it's a subsitution error
            if (queryRow[i] != refRow[i])
            {
                matchCount = 0;
                continue;
            }
            
            // If we made it this far it's a match - increment the counter
            matchCount++;

            // If we've seen enough matches in a row, return our position
            if (matchCount >= minMatchLength_)
            {
                return i + minMatchLength_;
            }
        }

        // If we made it this far, something went wrong - return an error
        return -1;
    }

    /// Since our alignments are local (i.e. relative to the position of the
    /// infixes used) we need to calculate and store their global position
    /// in a useful way.  This function recycles SeqAn's build in SimpleSeed
    /// class to do just that.
    ///
    /// \param  region  The TSeed object to store the alignment location in
    /// \param  align  The alignment object for which to record it's location
    /// \param  queryAlignStart  The number of bases trimmed from the infix
    ///                          that became the first alignment row
    /// \param  refAlignStart  The number of bases trimmed from the infix
    ///                        that became the second alignment row
    ///
    void SetAlignmentRegion(const TAlign& align, 
                            const size_t queryAlignStart, 
                            const size_t refAlignStart,
                            TSeed* region)
    {
        using namespace seqan;
        setBeginPositionH(*region, queryAlignStart + beginPosition(row(align, 0)));
        setEndPositionH(*region,   queryAlignStart + endPosition(row(align, 0)));
        setBeginPositionV(*region, refAlignStart + beginPosition(row(align, 1)));
        setEndPositionV(*region,   refAlignStart + endPosition(row(align, 1)));
    }

    /// We wish to only consider alignments that cover at least a certain
    /// span in the query.
    /// 
    /// \param  row  The sequence in which we wish to know the alignment span
    /// \param  alignStart  The start position, with respect to the alignment
    /// \param  alignEnd  The end position, with respect to the alignment
    //
    /// \return  size_t  The distance between the alignment start and end,
    ///                  with respect to the sequence provided
    inline
    size_t LengthInSequence(const TRow& row,
                            const int alignStart,
                            const int alignEnd)
    {
        using namespace seqan;

        return toSourcePosition(row, alignEnd) - toSourcePosition(row, alignStart);
    }

    /// Is a given character counted as a gap or not?
    ///
    /// \param  base  The character in question
    /// \return  bool  Whether the input is a gap character
    inline 
    bool IsGap(const char base)
    {
        return base == gapValue_;
    }

private:  // User-defined member variables
    const TStringSet& references_;
    const size_t flankingBuffer_;
    const size_t minQueryLength_;
    const size_t minMatchLength_;
    const size_t bandExtension_;
   
private:  // Class-defined memeber variables
    static constexpr float minAccuracyImprovement_ = 3.0f;
    
    // These aren't static since C doesn't allow static function calls
    const TScoringScheme scoringScheme_ = TScoringScheme(5, -6, -4, -4, -5, -5);
    const TAlignConfig alignConfig_ = TAlignConfig();
    const char gapValue_ = seqan::gapValue<char>();

};

} // SparseAlignment
} // PacBio
