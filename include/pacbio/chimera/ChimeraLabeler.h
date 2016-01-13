// Copyright (c) 2015, Pacific Biosciences of California, Inc.
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

// Author: Brett Bowman

#pragma once

#include <limits>
#include <vector>
#include <set>
#include <iostream>
#include <algorithm>
#include <string>

#include <seqan/sequence.h>
#include <seqan/graph_msa.h>

#include "ChimeraLabel.h"

namespace PacBio {
namespace Chimera {
    
namespace {

    /// Seprates a string on a specified delimiter
    ///
    /// \param s      Input string
    /// \param delim  Delimiter character
    ///
    /// \return Vector of sub-strings of the input string
    static std::vector<std::string> Split(const std::string& s, char delim)
    {
        std::vector<std::string> elems;
        std::stringstream ss(s);
        std::string item;
        while (std::getline(ss, item, delim)) {
            elems.push_back(item);
        }
        return elems;
    }

} // Anonymous namespace

///
/// Chimera detector - this is an implementation of the UCHIME algorithm, with a few generalizations:
/// Edgar, Robert C., et al. "UCHIME improves sensitivity and speed of chimera detection." Bioinformatics 27.16 (2011): 2194-2200.
///
/// Notes: We do full length alignments between the database and the test sequence, then find the best scoring splice site for each pair of database - test seq alignments
/// This is probably not as scalable as doing chunkwise alignments and only considering the best alignment to each chunk, but it's presumably more sensitive...
/// May need to speed this up (and lot's of other things) if the number of haplotypes gets very large.
///
class ChimeraLabeler
{
public:  // structors
    // Default constructor
    explicit ChimeraLabeler(double minChimeraScoreArg,
                            bool verbose)
            : minChimeraScore(minChimeraScoreArg)
            , verbose_(verbose)
            {};
    // Move constructor
    ChimeraLabeler(ChimeraLabeler&& src) = delete;
    // Copy constructor
    ChimeraLabeler(const ChimeraLabeler& src) = delete;
    // Move assignment constructor
    ChimeraLabeler& operator=(ChimeraLabeler&& rhs) = delete;
    // Copy assignment constructor
    ChimeraLabeler& operator=(const ChimeraLabeler& rhs) = delete;
    // Destructor
    ~ChimeraLabeler() = default;

private:  // Type definitions
    typedef seqan::Align<seqan::Dna5String, seqan::ArrayGaps> TAlign;
    typedef seqan::Score<int32_t, seqan::Simple> TScore;
    typedef seqan::Row<TAlign>::Type TRow;
    typedef seqan::Iterator<TRow>::Type TRowIter;
    typedef seqan::IterMakeConst<TRowIter>::Type TRowIterC;
    enum orientation_t { A_B, B_A, NA };

private:  // Instance variables
    const double minChimeraScore;
    const uint32_t beta = 4;
    const double pseudocount = 2.0f;
    const uint32_t chunks = 4;
    const bool verbose_;
    const TScore scoringScheme_ = TScore(2, -5, -3, -3);

private:  // State variables
    std::vector<std::string> ids_;
    std::vector<seqan::Dna5String> nonChimeras_;
    uint32_t minSize_   = std::numeric_limits<uint32_t>::max();
    size_t numAnalyzed_ = 0;

public:  // Modifying methods

    void clear()
    {
        ids_.clear();
        nonChimeras_.clear();
        minSize_     = std::numeric_limits<uint32_t>::max();
        numAnalyzed_ = 0;
    }

    /// \brief Label a vector of sequence records as Chimeric or not.
    ///        Secondary entry-point.
    ///
    /// \param A vector of all of the available sequence ids as strings
    /// \param A vector of all of the available sequences as Dna5Strings
    ///
    /// \return A set of labels representing the chimeric parents (if any) for
    ///         each input sequence
    ///
    std::vector<ChimeraLabel> LabelChimeras(const std::vector<std::string>& idList,
                                            const std::vector<std::string>& seqList,
                                            const std::vector<uint32_t>& sizeList)
    {
        std::vector<seqan::Dna5String> dnaStringList;
        for (const auto& seq : seqList)
            dnaStringList.emplace_back(seq);
        return LabelChimeras(idList, dnaStringList, sizeList);
    }

    /// \brief Label a vector of sequence records as Chimeric or not.
    ///        Secondary entry-point.
    ///
    /// \param A vector of all of the available sequence ids as strings
    /// \param A vector of all of the available sequences as Dna5Strings
    ///
    /// \return A set of labels representing the chimeric parents (if any) for
    ///         each input sequence
    ///
    std::vector<ChimeraLabel> LabelChimeras(const std::vector<std::string>& idList,
                                            const std::vector<std::string>& seqList)
    {
        std::vector<uint32_t> sizeList = ParseNumReads(idList);
        std::vector<seqan::Dna5String> dnaStringList;
        for (const auto& seq : seqList)
            dnaStringList.emplace_back(seq);
        return LabelChimeras(idList, dnaStringList, sizeList);
    }

    /// \brief Label a vector of sequence records as Chimeric or not.
    ///        Secondary entry-point.
    ///
    /// \param A vector of all of the available sequence ids as strings
    /// \param A vector of all of the available sequences as Dna5Strings
    ///
    /// \return A set of labels representing the chimeric parents (if any) for
    ///         each input sequence
    ///
    std::vector<ChimeraLabel> LabelChimeras(const std::vector<std::string>& idList,
                                            const std::vector<seqan::Dna5String>& seqList)
    {
        std::vector<uint32_t> sizeList = ParseNumReads(idList);
        return LabelChimeras(idList, seqList, sizeList);
    }

    /// \brief Label a vector of sequence records as Chimeric or not.
    ///        Main entry-point.
    ///
    /// \param A vector of all of the available sequence ids as strings
    /// \param A vector of all of the available sequences as Dna5Strings
    /// \param A vector of all of the available sequence "sizes" or prevalence levels
    ///
    /// \return A set of labels representing the chimeric parents (if any) for
    ///         each input sequence
    ///
    std::vector<ChimeraLabel> LabelChimeras(const std::vector<std::string>& ids,
                                            const std::vector<seqan::Dna5String>& seqs,
                                            const std::vector<uint32_t>& sizes)
    {
        // Declare the output vector now
        std::vector<ChimeraLabel> output;

        // Iterate over each Fasta record in order of their size
        for (uint32_t i = 0; i < ids.size(); ++i)
        {
            ChimeraLabel label = LabelChimera(ids[i], seqs[i], sizes[i]);
            output.push_back(label);
        }

        return output;
    }

    ChimeraLabel LabelChimera(const std::string& id,
                              const seqan::Dna5String& sequence,
                              const uint32_t size)
    {
        ChimeraLabel label(id);

        // Error-out if sequences are presented out of order
        if (size > minSize_)
            throw std::runtime_error("Sequences analyzed out of order!");

        // Declare containers for tracking non-Chimeric parents
        std::vector<uint32_t> parentIds;

        // First two sequences do not have enough parents, assumed real
        if (ids_.size() < 2)
        {
            if (verbose_)
                std::cout << "Consensus '" << id << "' is abundant, assumed real" << std::endl;
#ifdef PBLOG_INFO
            PBLOG_INFO << "Consensus '" << id << "' is abundant, assumed real";
#endif
            // Create a default label for the assumed-non-chimeric read
            AddNonChimera(id, sequence, size);
            return label;
        }
        else  // Otherwise we align
        {
            // Find probable parents from the highest scoring SW templates
            parentIds = FindParents(sequence);

            // If there is only one high scoring parent, then the sequence
            //    probably represents a true allele and we keep it.
            if (parentIds.size() == 1)
            {
                if (verbose_)
                    std::cout << "Consensus '" << id << "' has only one proposed parent, assumed real" << std::endl;
#ifdef PBLOG_INFO
                PBLOG_INFO << "Consensus '" << id << "' has only one proposed parent, assumed real";
#endif
                // Add a default label for the non-chimeric read
                AddNonChimera(id, sequence, size);
                return label;
            }
            else  // Otherwise we need to test it for chimerism
            {
                auto label = TestPossibleChimera(id, sequence, parentIds);
                
                if (verbose_)
                    std::cout << "Consensus '" << id << "' has a possible cross-over at " 
                              << label.crossover << " with a score of " << label.score << std::endl;
                    std::cout << "Possible parents are '" << label.leftParentId << "' and '"
                              << label.rightParentId << "'" << std::endl;
#ifdef PBLOG_INFO
                    PBLOG_INFO << "Consensus '" << id << "' has a possible cross-over at " 
                               << label.crossover << " with a score of " << label.score;
                    PBLOG_INFO << "Possible parents are '" << label.leftParentId << "' and '"
                               << label.rightParentId << "'";
#endif

                // If the score is high enough, set the flag 
                if (label.score > minChimeraScore) {
                    label.chimeraFlag = true;
                } else {
                    AddNonChimera(id, sequence, size);
                }

                return label;
            }
        }
    }


private:  // non-modifying methods

    void AddNonChimera(const std::string& id,
                       const seqan::Dna5String sequence,
                       const uint32_t size)
    {
        ids_.push_back(id);
        nonChimeras_.push_back(sequence);
        AddReverseComplement(sequence);
        minSize_ = std::min(size, minSize_);
    }

    void AddReverseComplement(const seqan::Dna5String& sequence)
    {
        seqan::Dna5String revComSeq = sequence;
        seqan::reverseComplement(revComSeq);
        nonChimeras_.push_back(revComSeq);
    }

    /// \brief Append all of the reverse-complements to the end of a vector of sequences
    ///
    /// \param A vector of N DNA sequences
    ///
    /// \return A vector of 2*N DNA sequences
    ///
    void AddReverseComplements(std::vector<seqan::Dna5String>& seqList)
    {
        for (const auto& sequence : seqList)
        {
            auto revComSeq = sequence;
            seqan::reverseComplement(revComSeq);
            seqList.push_back(revComSeq);
        }
        return;
    }

    /// \brief Find the most probable parents for a possible chimera
    /// 
    /// \param A vector of all of the available sequences
    /// \param A vector of indices of non-chimeric sequences in that vector
    /// \param The index of the current sequence to be tested
    ///
    /// \return A set of indices representing the best scoring parents for the
    ///         various regions of the current sequences
    ///
    std::vector<uint32_t> FindParents(const seqan::Dna5String& sequence)
    {
        // Declare the output variable and the set we will build it from
        std::vector<uint32_t> output;
        std::set<uint32_t> parentIds;

        // Initialize the pairwise-alignment object we will be reusing
        TAlign align;
        seqan::resize(rows(align), 2);

        // Pre-calculate the size of each chunk
        uint32_t chunkSize = seqan::length(sequence) / chunks;

        // Iterate over each chunk, aligning it to all possible parents
        for (uint32_t i = 0; i < chunks; ++i)
        {
            // Initialize the alignment with the current sequence chunk
            uint32_t chunkStart = i * chunkSize;
            uint32_t chunkEnd   = chunkStart + chunkSize;
            const auto chunkSeq = infix(sequence, chunkStart, chunkEnd);
            seqan::assignSource(seqan::row(align, 0), chunkSeq);

            // Initialize loop variables;
            uint32_t score      = 0;
            uint32_t maxScore   = 0;
            uint32_t maxParent  = 0;

            // iterate over each non-Chimeric sequence
            for (size_t i = 0; i < nonChimeras_.size(); ++i)
            {
                // Fill out the alignment with the current parents
                seqan::assignSource(seqan::row(align, 1), nonChimeras_[i]);
                score = seqan::localAlignment(align, scoringScheme_);

                // If the current parent is better than the best, keep it
                if (score > maxScore)
                {
                    maxScore = score;
                    maxParent = i;
                }
            }

            // Add the best parent for this chunk to the set
            parentIds.insert(maxParent);
        }

        // Convert the set of parents to a vector for down-stream use
        std::move(parentIds.begin(), parentIds.end(), std::back_inserter(output));
        return output;
    }

    /// \brief Identify the highest-scoring chimeric explaination for a query
    ///        from a list of
    ///
    /// \param A pointer to a vector of sequences
    /// \param A pointer to a vector of sequence ids
    /// \param An unsigned int for the index of the query in the above vectors
    /// \param A vector of unsigned ints for
    ///
    /// \return A ChimeraLabel for the highest-scoring chimeric explaination
    ///         for a given query
    ///
    ChimeraLabel TestPossibleChimera(
            const std::string& id,
            const seqan::Dna5String& sequence,
            const std::vector<uint32_t>& possibleParents)
    {
        // First we align the query to all of the possible parents
        //  auto alignments = GetFullAlignments(sequences, index, 
        //     possibleParents);
        auto alignments = GetMultiSequenceAlignment(sequence, possibleParents);

        // Initialize two Labels - one for the max and one for the loop variable
        ChimeraLabel bestLabel(id);

        // Loop variables for the names of the identified parents
        const uint32_t idCount = ids_.size();
        std::string parentA;
        std::string parentB;

        // Iterate over all possible "Parent A"'s
        for (uint32_t i = 1; i < possibleParents.size(); ++i)
        {
            uint32_t parentAIdx = possibleParents[i];
            parentA = ids_[parentAIdx/2];

            // Iterate over all possible "Parent B"'s
            for (uint32_t j = 0; j < i; ++j)
            {
                uint32_t parentBIdx = possibleParents[j];
                parentB = ids_[parentBIdx/2];

                // For a given ParentA and ParentB, what is the maximum
                //     possible Chimera score?
                ChimeraLabel label = ScorePossibleChimera(alignments, id,
                                                          parentA, parentB,
                                                          i+1, j+1);

                // Keep the highest scoring label
                if (label.score > bestLabel.score)
                    bestLabel = std::move(label);
            }
        }

        return bestLabel;
    }

    /// \brief Generate an MSA of a query sequence and all possible parents
    ///
    /// \param A pointer to a vector of possible sequences
    /// \param An unsigned int for the index of the query in the vector of sequences
    /// \param A vector of unsigned ints for the indices of all possible parents
    ///
    /// \return a pointer to an MSA
    ///
    TAlign GetMultiSequenceAlignment(
            const seqan::Dna5String& sequence,
            const std::vector<uint32_t>& parentIds)
    {
        // Initialize the alignment with the query sequence
        TAlign align;
        seqan::resize(rows(align), parentIds.size()+1);
        seqan::assignSource(seqan::row(align, 0), sequence);

        // Successively add each parent to the alignment
        for (uint32_t i = 1; i <= parentIds.size(); ++i)
        {
            uint32_t parentIdx = parentIds[i-1];
            seqan::assignSource(seqan::row(align, i), nonChimeras_[parentIdx]);
        }

        // Perform the alignment operation and return
        seqan::globalMsaAlignment(align, scoringScheme_);
        return align;
    }

    /// \brief Scan an MSA of sequences for all possible chimeric break-points
    ///        that could explain the query as a composite of the parents
    ///
    /// \param A pointer to an MSA of sequences
    /// \param A string for the name of the Query
    /// \param A string for the name of the first parent
    /// \param A string for the name of the second parent
    /// \param An unsigned int for the row-index of the first parent
    /// \param An unsigned int for the row-index of the second parent
    ///
    /// \return A ChimeraLabel of the highest-scoring possible chimera
    ///
    ChimeraLabel ScorePossibleChimera(
            const TAlign& alignment,
            const std::string& queryId,
            const std::string& parentAId,
            const std::string& parentBId,
            const uint32_t firstIdx,
            const uint32_t secondIdx)
    {
        // Second, extract references to the rows we need to inspect
        const TRow& queryRow = seqan::row(alignment,         0);
        const TRow& parentA  = seqan::row(alignment,  firstIdx);
        const TRow& parentB  = seqan::row(alignment, secondIdx);

        // Initialize const-iterators for each row
        TRowIterC itQ = seqan::begin(queryRow);
        TRowIterC itA = seqan::begin(parentA);
        TRowIterC itB = seqan::begin(parentB);

        // Initialize end-point for the const-iterators (since we're traversing
        //    an MSA all 3 iterators will terminate at the same point)
        TRowIterC qEnd = seqan::end(queryRow);

        // Count variables
        uint32_t rightA   = 0;
        uint32_t rightB   = 0;
        uint32_t rightAbs = 0;
        uint32_t leftA    = 0;
        uint32_t leftB    = 0;
        uint32_t leftAbs  = 0;

        // First we iterate once, counting up the total number of A/B/Abs votes
        //    to initialize the counts for the right-hand segment
        for (; itQ != qEnd; ++itQ, ++itA, ++itB)
        {
            // We can't trust our consensus around gaps for low-coverage reads
            //     so we skip them
            if (seqan::isGap(itQ) || seqan::isGap(itA) || seqan::isGap(itB))
                continue;

            // Identical bases are useless for Chimera Detection - skip
            if (*itQ == *itA && *itQ == *itB)
                continue;

            // If we made it here, count the difference as a vote
            if (*itA == *itB) {
                ++rightAbs;
            } else if (*itQ == *itA) {
                ++rightA;
            } else if (*itQ == *itB) {
                ++rightB;
            }
        }

        // Re-initialize iterators for each row
        itQ = seqan::begin(queryRow);
        itA = seqan::begin(parentA);
        itB = seqan::begin(parentB);

        // Initialize variables for the maximum chimera
        double maxChimeraScore = 0.0f;
        uint32_t maxChimeraCrossover = 0;
        orientation_t maxChimeraOrientation = NA;

        // And per-iteration variables
        double abScore, baScore, chimeraScore;
        uint32_t chimeraCrossover = 0;
        orientation_t chimeraOrientation;

        // Second time we iterate, we subtract from the right
        for (; itQ != qEnd; ++itQ, ++itA, ++itB)
        {
            // If the Query sequence is not at a gap, increment the
            //     tracking variable for the crossover postion
            if (seqan::isGap(itQ))
                continue;
            else
                ++chimeraCrossover;

            // We can't trust our consensus around gaps for low-coverage reads
            //     so we skip them
            if (seqan::isGap(itA) || seqan::isGap(itB))
                continue;

            // Identical bases are useless for Chimera Detection - skip
            if (*itQ == *itA && *itQ == *itB)
                continue;

            // If we made it here, count the differences as a vote
            if (*itA == *itB) {
                ++leftAbs;
                --rightAbs;
            } else if (*itQ == *itA) {
                ++leftA;
                --rightA;
            } else if (*itQ == *itB) {
                ++leftB;
                --rightB;
            }

            // If we've exhausted the right-side votes, there are no more
            //     possible break points
            if (rightA == 0 && rightB == 0)
                break;

            // If haven't seen any left-side votes, we haven't reached any
            //     possible break points yet
            if (leftA == 0 && leftB == 0)
                continue;

            // If the Left leans A and the right leans B, test "AAABBB"
            if (leftA > leftB && rightA < rightB) 
            {
                chimeraScore = ScoreBreakPoint(leftA, leftB, leftAbs,
                        rightB, rightA, rightAbs);
                chimeraOrientation = A_B;
            }  // If the Left leans B and the right leans A, test "BBBAAA"
            else if (leftA < leftB && rightA > rightB) 
            {
                chimeraScore = ScoreBreakPoint(leftB, leftA, leftAbs,
                        rightA, rightB, rightAbs);
                chimeraOrientation = B_A;
            }  // If either left or right isn't clear, test both options
            else 
            {
                abScore = ScoreBreakPoint(leftA, leftB, leftAbs,
                        rightB, rightA, rightAbs);
                baScore = ScoreBreakPoint(leftB, leftA, leftAbs,
                        rightA, rightB, rightAbs);
                if (abScore > baScore) {
                    chimeraScore = abScore;
                    chimeraOrientation = A_B;
                } else {
                    chimeraScore = baScore;
                    chimeraOrientation = B_A;
                }
            }

            // Keep the best chimera we've seen so far
            if (chimeraScore > maxChimeraScore)
            {
                maxChimeraScore = chimeraScore;
                maxChimeraCrossover = chimeraCrossover;
                maxChimeraOrientation = chimeraOrientation;
            }
        }

        // Build and return an appropriately oriented ChimeraLabel
        if (maxChimeraOrientation == A_B)
            return ChimeraLabel(queryId, parentAId, parentBId,
                    maxChimeraCrossover, maxChimeraScore);
        else
            return ChimeraLabel(queryId, parentBId, parentAId,
                    maxChimeraCrossover, maxChimeraScore);
    }

public:  // Public Static Methods

    static uint32_t ParseNumReads(const std::string id)
    {
        const auto& parts = Split(id, '_');
        const auto& numReadsString = parts[3].substr(8);
        const uint32_t numReads = std::stoi(numReadsString);
        return numReads;
    }

    static std::vector<uint32_t> ParseNumReads(const std::vector<std::string> ids)
    {
        std::vector<uint32_t> retval;
        for (size_t i = 0; i < ids.size(); ++i)
            retval.emplace_back(ParseNumReads(ids[i]));
        return retval;
    }

private:  // Inlined methods

    /// \brief Calculates the H-score for a chimeric alignment as per Edgar(2011)
    ///
    /// \param An unsigned int of votes for left-side similarity
    /// \param An unsigned int of votes against left-side similarity
    /// \param An unsigned int of votes neither for or against the left side
    /// \param An unsigned int of votes for right-side similarity
    /// \param An unsigned int of votes against right-side similarity
    /// \param An unsigned int of votes neither for or against the right side
    ///
    /// \return A double representing the strength of the similarity between the query
    ///         sequence and a hypothetical chimera composed of left and right segments
    ///         taken from existing sequences
    ///
    inline double ScoreBreakPoint(const uint32_t leftYesVotes,
                                  const uint32_t leftNoVotes,
                                  const uint32_t leftAbsVotes,
                                  const uint32_t rightYesVotes,
                                  const uint32_t rightNoVotes,
                                  const uint32_t rightAbsVotes)
    {
        // Score the left and right segments independently
        double leftScore = ScoreSegment(leftYesVotes, leftNoVotes,
                leftAbsVotes);
        double rightScore = ScoreSegment(rightYesVotes, rightNoVotes,
                rightAbsVotes);

        // Return their product
        return leftScore * rightScore;
    }

    /// \brief Calculates the H-score for a pairwise alignment segment,
    ///        as per Edgar(2011)
    ///
    /// \param An unsigned int of votes for similarity
    /// \param An unsigned int of votes against similarity
    /// \param An unsigned int of votes neither for nor against
    ///
    /// \return A double representing the strength of the similarity between two
    ///         sequences
    ///
    inline double ScoreSegment(const uint32_t yesVotes,
                               const uint32_t noVotes,
                               const uint32_t absVotes)
    {
        return yesVotes / (beta * (noVotes + pseudocount) + absVotes);
    }
};

}  // namespace Chimera
}  // namespace PacBio
