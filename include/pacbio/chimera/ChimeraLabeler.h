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

#include <seqan/sequence.h>
#include <seqan/graph_msa.h>

#include <vector>
#include <set>
#include <iostream>
#include <algorithm>
#include <string>

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
    explicit ChimeraLabeler(double minChimeraScoreArg)
            : minChimeraScore(minChimeraScoreArg)
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
    const double minParentScore = 1.0f;
    const uint32_t beta = 4;
    const double pseudocount = 2.0f;
    const uint32_t chunks = 4;
    const TScore scoringScheme = TScore(2, -5, -3, -3);

public:  // non-modifying methods

    /// \brief Label a vector of sequence records as Chimeric or not.
    ///        Secondary entry-point.
    ///
    /// \param A vector of all of the available sequence ids as strings
    /// \param A vector of all of the available sequences as Dna5Strings
    ///
    /// \return A set of labels representing the chimeric parents (if any) for
    ///         each input sequence
    ///
    std::vector<ChimeraLabel> Label(const std::vector<std::string>& idList,
                                    const std::vector<std::string>& seqList)
    {
        std::vector<uint32_t> sizeList = ParseNumReads(idList);
        std::vector<seqan::Dna5String> dnaStringList;
        for (const auto& seq : seqList)
            dnaStringList.emplace_back(seq);
        return Label(idList, dnaStringList, sizeList);
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
    std::vector<ChimeraLabel> Label(const std::vector<std::string>& idList,
                                    const std::vector<seqan::Dna5String>& seqList)
    {
        std::vector<uint32_t> sizeList = ParseNumReads(idList);
        return Label(idList, seqList, sizeList);
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
    std::vector<ChimeraLabel> Label(const std::vector<std::string>& idList,
                                    const std::vector<seqan::Dna5String>& seqList2,
                                    const std::vector<uint32_t>& sizeList)
    {
        // Declare the output vector now
        std::vector<ChimeraLabel> output;

        // Sanity-check the sizes of the vectors we've been give.  All three should be
        // of the same size
        if (idList.size() != seqList2.size())
            throw std::runtime_error("Sizes of Id List and Sequence List do not match!");
        if (idList.size() != sizeList.size())
            throw std::runtime_error("Sizes of Id List and Size List do not match!");

        // Sanity-check the ordering of the sizes.
        for (size_t i = 1; i < sizeList.size(); ++i)
            if (sizeList[i] > sizeList[i-1])
                throw std::runtime_error("Sequences are not sorted by size!");

        // Add the reverse complements for each sequence
        std::vector<seqan::Dna5String> seqList;
        for (const auto& sequence : seqList2)
            seqList.push_back(sequence);
        AddReverseComplements(seqList);

        // Declare loop variables
        uint32_t N = idList.size();
        std::string id;
        seqan::Dna5String sequence;

        // Declare containers for tracking non-Chimeric parents
        std::vector<uint32_t> nonChimeras;
        std::vector<uint32_t> parentIds;

        // Iterate over each Fasta record in order of their size
        for (uint32_t i = 0; i < idList.size(); ++i)
        {
            // Pull out the current Id / Sequence / Size
            id       = idList[i];
            sequence = seqList[i];

            std::cerr << "Analyzing sequence #" << i + 1
                      << " of " << idList.size()
                      << ", '" << id << "', supported by "
                      << sizeList[i] << " reads"
                      << std::endl;

            // First two sequences do not have enough parents, assumed real
            if (output.size() < 2)
            {
                // Create a default label for the assumed-non-chimeric read
                output.emplace_back(id);
            }
            else  // Otherwise we align
            {
                // Find probable parents from the highest scoring SW templates
                parentIds = FindParents(seqList, nonChimeras, i);

                // If there is only one high scoring parent, then the sequence
                //    probably represents a true allele and we keep it.
                if (parentIds.size() == 1)
                {
                    // Add a default label for the non-chimeric read
                    output.emplace_back(id);
                }
                else  // Otherwise we need to test it for chimerism
                {
                    auto label = TestPossibleChimera(seqList, idList, i, 
                            parentIds);

                    // If the score is high enough, set the flag and save it
                    if (label.score > minChimeraScore) {
                        label.chimeraFlag = true;
                        output.push_back(label);
                    } else { // Otherwise add the default label
                        output.emplace_back(id);
                    }
                }
            }

            // Add the current index and it's RC to known non-chimeras
            nonChimeras.push_back(i);
            nonChimeras.push_back(i+N);
        }

        return output;  // Implicit move semantics take care that ChimeraLabels are moved.
    }

private:
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
    std::vector<uint32_t> FindParents(const std::vector<seqan::Dna5String>& sequences,
                                      const std::vector<uint32_t>& nonChimericIdx,
                                      const uint32_t index)
    {
        // Declare the output variable and the set we will build it from
        std::vector<uint32_t> output;
        std::set<uint32_t> parentIds;

        // Initialize the pairwise-alignment object we will be reusing
        TAlign align;
        seqan::resize(rows(align), 2);

        // Pre-calculate the size of each chunk
        uint32_t chunkSize = seqan::length(sequences[index]) / chunks;

        // Iterate over each chunk, aligning it to all possible parents
        for (uint32_t i = 0; i < chunks; ++i)
        {
            // Initialize the alignment with the current sequence chunk
            uint32_t chunkStart = i * chunkSize;
            uint32_t chunkEnd   = chunkStart + chunkSize;
            const auto chunkSeq = infix(sequences[index], chunkStart, chunkEnd);
            seqan::assignSource(seqan::row(align, 0), chunkSeq);

            // Initialize loop variables;
            uint32_t score      = 0;
            uint32_t maxScore   = 0;
            uint32_t testParent = 0;
            uint32_t maxParent  = 0;

            // iterate over each non-Chimeric sequence
            for (size_t i = 0; i < nonChimericIdx.size(); ++i)
            {
                // Fill out the alignment with the current parents
                testParent = nonChimericIdx[i];
                seqan::assignSource(seqan::row(align, 1), 
                        sequences[testParent]);
                score = seqan::localAlignment(align, scoringScheme);

                // If the current parent is better than the best, keep it
                if (score > maxScore)
                {
                    maxScore = score;
                    maxParent = testParent;
                }
            }

            // Add the best parent for this chunk to the set
            parentIds.insert(maxParent);
        }

        // Convert the set of parents to a vector for down-stream use
        std::move(parentIds.begin(), parentIds.end(), 
                  std::back_inserter(output));
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
            const std::vector<seqan::Dna5String>& sequences,
            const std::vector<std::string>& ids,
            const uint32_t index,
            const std::vector<uint32_t>& possibleParents)
    {
        // First we align the query to all of the possible parents
        //  auto alignments = GetFullAlignments(sequences, index, 
        //     possibleParents);
        auto alignments = GetMultiSequenceAlignment(sequences, index, 
                possibleParents);

        // Initialize two Labels - one for the max and one for the loop variable
        std::string id = ids[index];
        ChimeraLabel bestLabel(id);
        ChimeraLabel label;

        // Loop variables for the names of the identified parents
        const uint32_t idCount = ids.size();
        std::string parentA;
        std::string parentB;

        // Iterate over all possible "Parent A"'s
        for (uint32_t i = 1; i < possibleParents.size(); ++i)
        {
            uint32_t parentAIdx = possibleParents[i];
            if (parentAIdx >= idCount)
                parentA = ids[parentAIdx-idCount];
            else
                parentA = ids[parentAIdx];

            // Iterate over all possible "Parent B"'s
            for (uint32_t j = 0; j < i; ++j)
            {
                uint32_t parentBIdx = possibleParents[j];
                if (parentBIdx >= idCount)
                    parentB = ids[parentBIdx-idCount];
                else
                    parentB = ids[parentBIdx];

                // For a given ParentA and ParentB, what is the maximum
                //     possible Chimera score?
                label = ScorePossibleChimera(alignments, id,
                                             parentA, parentB,
                                             i+1, j+1);

                // Keep the highest scoring label
                if (label.score > bestLabel.score)
                {
                    bestLabel = std::move(label);
                }
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
            const std::vector<seqan::Dna5String>& sequences,
            const uint32_t index,
            const std::vector<uint32_t>& parentIds)
    {
        // Initialize the alignment with the query sequence
        TAlign align;
        seqan::resize(rows(align), parentIds.size()+1);
        seqan::assignSource(seqan::row(align, 0), sequences[index]);

        // Successively add each parent to the alignment
        for (uint32_t i = 1; i <= parentIds.size(); ++i)
        {
            uint32_t parentIdx = parentIds[i-1];
            seqan::assignSource(seqan::row(align, i), 
                    sequences[parentIdx]);
        }

        // Perform the alignment operation and return
        seqan::globalMsaAlignment(align, scoringScheme);
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
