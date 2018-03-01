// Author: Brett Bowman

#pragma once

#include <cassert>

#include <algorithm>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <pbcopper/utility/StringUtils.h>

#include <pbbam/Cigar.h>

#include <pacbio/align/LocalAlignment.h>
#include <pacbio/data/Sequence.h>

#include "ChimeraLabel.h"

namespace PacBio {
namespace Chimera {

///
/// Chimera detector - this is an implementation of the UCHIME algorithm, with a
/// few generalizations:
/// Edgar, Robert C., et al. "UCHIME improves sensitivity and speed of chimera
/// detection." Bioinformatics 27.16 (2011): 2194-2200.
///
/// Notes: We do full length alignments between the database and the test
/// sequence, then find the best scoring splice site for each pair of database -
/// test seq alignments
/// This is probably not as scalable as doing chunkwise alignments and only
/// considering the best alignment to each chunk, but it's presumably more
/// sensitive...
/// May need to speed this up (and lot's of other things) if the number of
/// haplotypes gets very large.
///
class ChimeraLabeler
{
public:  // structors
    // Default constructor
    explicit ChimeraLabeler(double minChimeraScoreArg = 1.0, size_t maxChimeraSupportArg = 100,
                            bool verboseArg = false)
        : minChimeraScore(minChimeraScoreArg)
        , maxChimeraSupport(maxChimeraSupportArg)
        , chunks(4)
        , verbose(verboseArg){};
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
    enum orientation_t
    {
        A_B,
        B_A,
        NA
    };

private:  // Instance variables
    const double minChimeraScore;
    const size_t maxChimeraSupport;
    const size_t chunks;
    const bool verbose;

private:  // State variables
    std::vector<std::string> ids_;
    std::vector<std::string> nonChimeras_;
    size_t minSize_ = std::numeric_limits<size_t>::max();
    size_t numAnalyzed_ = 0;

public:  // Modifying methods
    /// \brief Clear and reset all stateful variables to default
    ///
    void clear()
    {
        ids_.clear();
        nonChimeras_.clear();
        minSize_ = std::numeric_limits<size_t>::max();
        numAnalyzed_ = 0;
    }

    /// \brief Label a vector of sequence records as Chimeric or not.
    ///        Secondary entry-point.
    ///
    /// \param A vector of all of the sequence ids as strings
    /// \param A vector of all of the sequences as strings
    ///
    /// \return A set of labels representing the chimeric parents (if any) for
    ///         each input sequence
    ///
    std::vector<ChimeraLabel> LabelChimeras(const std::vector<std::string>& idList,
                                            const std::vector<std::string>& seqList)
    {
        const auto sizeList = ParseNumReads(idList);
        return LabelChimeras(idList, seqList, sizeList);
    }

    /// \brief Label a vector of sequence records as Chimeric or not.
    ///        Main entry-point.
    ///
    /// \param A vector of all the sequence ids as strings
    /// \param A vector of all the sequences as strings
    /// \param A vector of integers, representing the support for each sequence
    ///
    /// \return A set of labels representing the chimeric parents (if any) for
    ///         each input sequence
    ///
    std::vector<ChimeraLabel> LabelChimeras(const std::vector<std::string>& ids,
                                            const std::vector<std::string>& seqs,
                                            const std::vector<size_t>& sizes)
    {
        const size_t numSeqs = seqs.size();
        if (ids.size() != numSeqs || sizes.size() != numSeqs)
            throw std::runtime_error("Input containers must contain the same number of elements.");

        std::vector<ChimeraLabel> output;
        output.reserve(numSeqs);
        for (size_t i = 0; i < numSeqs; ++i) {
            output.push_back(LabelChimera(ids[i], seqs[i], sizes[i]));
        }
        return output;
    }

    /// \brief Label a vector of sequence records as Chimeric or not.
    ///        Main entry-point.
    ///
    /// \param The Id of the sequence as a string
    /// \param The sequence itself
    /// \param The support for the sequence, as an integer
    ///
    /// \return A label for whether and how the query is chimeric
    ///
    ChimeraLabel LabelChimera(const std::string& id, const std::string& sequence, const size_t size)
    {
        // Error-out if sequences are presented out of order
        if (size > minSize_) throw std::runtime_error("Sequences analyzed out of order!");

        // Declare containers for tracking non-Chimeric parents
        std::vector<size_t> parentIds;

        // The first two sequences do not have enough possible parents, and any
        //  records with a high enough number of supporting reads, assume real
        if (ids_.size() < 2 || size > maxChimeraSupport) {
            if (verbose)
                std::cout << "consensus '" << id << "' is abundant, assumed real" << std::endl;
#ifdef PBLOG_DEBUG
            PBLOG_DEBUG << "consensus '" << id << "' is abundant, assumed real";
#endif
            // Create a default label for the assumed-non-chimeric read
            AddNonChimera(id, sequence, size);
            return ChimeraLabel{id};
        } else  // Otherwise we align
        {
            // Find probable parents from the highest scoring SW templates
            parentIds = FindParents(sequence);

            // If there is only one high scoring parent, then the sequence
            //    probably represents a true allele and we keep it.
            if (parentIds.size() == 1) {
                if (verbose)
                    std::cout << "consensus '" << id
                              << "' has only one proposed parent, assumed real" << std::endl;
#ifdef PBLOG_DEBUG
                PBLOG_DEBUG << "consensus '" << id
                            << "' has only one proposed parent, assumed real";
#endif
                // Add a default label for the non-chimeric read
                AddNonChimera(id, sequence, size);
                return ChimeraLabel{id};
            } else  // Otherwise we need to test it for chimerism
            {
                auto label = TestPossibleChimera(id, sequence, parentIds);

                if (label.score > 0.0) {
                    if (verbose) {
                        std::cout << "consensus '" << id << "' has a possible cross-over at "
                                  << label.crossover << " with a score of " << label.score
                                  << std::endl;
                        std::cout << "possible parents are '" << label.leftParentId << "' and '"
                                  << label.rightParentId << "'" << std::endl;
                    }
#ifdef PBLOG_DEBUG
                    PBLOG_DEBUG << "consensus '" << id << "' has a possible cross-over at "
                                << label.crossover << " with a score of " << label.score;
                    PBLOG_DEBUG << "possible parents are '" << label.leftParentId << "' and '"
                                << label.rightParentId << "'";
#endif
                } else {
                    if (verbose)
                        std::cout
                            << "consensus '" << id
                            << "' had no possible cross-overs with positive scores, assumed real"
                            << std::endl;
#ifdef PBLOG_DEBUG
                    PBLOG_DEBUG
                        << "consensus '" << id
                        << "' has no possible cross-overs with positive scores, assumed real";
#endif
                }

                // If the score is high enough, set the flag, otherwise we add it to our
                // reference
                if (label.score > minChimeraScore) {
                    label.chimeraFlag = true;
                } else {
                    AddNonChimera(id, sequence, size);
                }

                return label;
            }
        }
    }

private:
    /// \brief Store a non-chimeric sequence to consider as a possible parent later
    ///
    /// \param The sequence Id as a string
    /// \param The sequence itself, as a string
    /// \param Its level of coverage or support
    ///
    void AddNonChimera(const std::string& id, const std::string& sequence, const size_t size)
    {
        ids_.push_back(id);
        nonChimeras_.push_back(sequence);
        nonChimeras_.push_back(PacBio::Data::ReverseComplement(sequence));
        minSize_ = std::min(size, minSize_);
    }

    /// \brief Find the most probable parents for a possible chimera
    ///
    /// \param the sequence to be tested
    ///
    /// \return A set of indices representing the best scoring parents for the
    ///         various regions of the current sequences
    ///
    std::vector<size_t> FindParents(const std::string& sequence) const
    {
        // Declare the output variable and the set we will build it from
        std::vector<size_t> output;
        std::set<size_t> parentIds;

        // Pre-calculate the size of each chunk
        const size_t chunkSize = sequence.size() / chunks;

        // Re-used alignment scoring scheme
        // TODO: Fix/Replace crashing SSW impl and revert/re-optimize the alignment params
        //const PacBio::Align::LocalAlignConfig alignConfig{2, 5, 3, 3};
        const PacBio::Align::LocalAlignConfig alignConfig =
            PacBio::Align::LocalAlignConfig::Default();

        // Iterate over each chunk, aligning it to all possible parents
        for (size_t i = 0; i < chunks; ++i) {
            // Initialize the alignment with the current sequence chunk
            const size_t chunkStart = i * chunkSize;
            const auto target = sequence.substr(chunkStart, chunkSize);

            // Initialize loop variables;
            size_t score = 0;
            size_t maxScore = 0;
            size_t maxParent = 0;

            // iterate over each non-Chimeric sequence
            for (size_t j = 0; j < nonChimeras_.size(); ++j) {
                // Fill out the alignment with the current parents
                const auto& query = nonChimeras_[j];

                // The underlying SSW impl finds a region in Seq2 to which to align Seq1.
                // This leads to banding problems and can crash if Seq2 is large or
                // repetitve, so we need to use the smaller sequence as Seq2 to avoid this.
                if (target.size() > query.size()) {
                    const auto al = PacBio::Align::LocalAlign(target, query, alignConfig);
                    score = al.Score();
                } else {
                    const auto al = PacBio::Align::LocalAlign(query, target, alignConfig);
                    score = al.Score();
                }

                // If the current parent is better than the best, keep it
                if (score > maxScore) {
                    maxScore = score;
                    maxParent = j;
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
    ///        from a list of possible parents
    ///
    /// \param The query Id as a string
    /// \param The query sequence as a string
    /// \param A vector of row indices indicating possible parent
    ///
    /// \return A ChimeraLabel for the highest-scoring chimeric explaination
    ///         for a given query
    ///
    ChimeraLabel TestPossibleChimera(const std::string& id, const std::string& sequence,
                                     const std::vector<size_t>& possibleParents) const
    {
        const auto alignments = GetMultiSequenceAlignment(sequence, possibleParents);

        // Initialize two Labels - one for the max and one for the loop variable
        ChimeraLabel bestLabel(id);

        // Loop variables for the names of the identified parents
        std::string parentA;
        std::string parentB;

        // Iterate over all possible "Parent A"'s
        for (size_t i = 1; i < possibleParents.size(); ++i) {
            size_t parentAIdx = possibleParents[i];
            parentA = ids_[parentAIdx / 2];

            // Iterate over all possible "Parent B"'s
            for (size_t j = 0; j < i; ++j) {
                size_t parentBIdx = possibleParents[j];
                parentB = ids_[parentBIdx / 2];

                // For a given ParentA and ParentB, what is the maximum
                //     possible Chimera score?
                ChimeraLabel label =
                    ScorePossibleChimera(alignments, id, parentA, parentB, i + 1, j + 1);

                // Keep the highest scoring label
                if (label.score > bestLabel.score && label.score > 0.0)
                    bestLabel = std::move(label);
            }
        }

        return bestLabel;
    }

    /// \brief Generate a (pseudo-)MSA of a query sequence and all possible parents
    ///
    /// \param target seequence
    /// \param A vector of unsigned ints for the indices of all possible parents
    ///
    /// \return a vector of aligned sequences
    ///
    std::vector<std::string> GetMultiSequenceAlignment(const std::string& targetSequence,
                                                       const std::vector<size_t>& parentIds) const
    {
        typedef std::map<size_t, uint32_t> InsertionMap;  // position -> numInsertions

        // store each parentId seq (both orientations) to align against target sequence
        std::vector<std::string> queries;
        queries.reserve(parentIds.size() * 2);
        for (const auto parentIdx : parentIds) {
            const auto& parentSeq = nonChimeras_[parentIdx];
            queries.emplace_back(parentSeq);
            queries.emplace_back(PacBio::Data::ReverseComplement(parentSeq));
        }

        // align all possible parents against target sequence
        const PacBio::Align::LocalAlignConfig alignConfig{2, 5, 3, 3};
        auto alignments = PacBio::Align::LocalAlign(targetSequence, queries, alignConfig);
        assert(alignments.size() == queries.size());
        assert(alignments.size() % 2 == 0);

        // setup alignment matrix, with gapped sequences, storing max insertions at each position
        std::vector<std::string> alignmentMatrix;
        alignmentMatrix.reserve(parentIds.size() + 1);
        alignmentMatrix.push_back(targetSequence);

        InsertionMap maxInsertionsLookup;
        std::vector<InsertionMap> perAlignmentInsertions;
        perAlignmentInsertions.push_back(InsertionMap{});  // empty map for target seq

        const auto numAlignments = alignments.size();
        for (size_t i = 0; i < numAlignments; i += 2) {
            assert((numAlignments - i) >= 2);

            const bool usingForward = (alignments.at(i).Score() >= alignments.at(i + 1).Score());
            const auto& align = (usingForward ? alignments.at(i) : alignments.at(i + 1));
            const auto& query = (usingForward ? queries.at(i) : queries.at(i + 1));
            const size_t alignBegin = static_cast<size_t>(align.TargetBegin());

            std::string gappedSeq;
            gappedSeq.reserve(targetSequence.size());

            // add gaps before alignment begin
            for (size_t j = 0; j < alignBegin; ++j)
                gappedSeq.push_back('-');

            size_t qPos = 0;
            size_t tPos = alignBegin;
            InsertionMap insertions;

            // fill out sequence, adding gaps for deletions
            const auto cigar = PacBio::BAM::Cigar::FromStdString(align.CigarString());
            for (const auto& op : cigar) {
                const auto type = op.Type();
                const auto length = op.Length();

                using PacBio::BAM::CigarOperationType;

                // store insertions for later
                if (type == CigarOperationType::INSERTION ||
                    type == CigarOperationType::SOFT_CLIP) {
                    const auto iter = maxInsertionsLookup.find(tPos);
                    if (iter == maxInsertionsLookup.cend() || (length > iter->second))
                        maxInsertionsLookup[tPos] = length;
                    insertions[tPos] = length;
                }

                for (size_t j = 0; j < length; ++j) {
                    switch (type) {
                        case CigarOperationType::SEQUENCE_MATCH:
                            assert(targetSequence.at(tPos) == query.at(qPos));  // fall-through
                        case CigarOperationType::ALIGNMENT_MATCH:               // fall-through
                        case CigarOperationType::SEQUENCE_MISMATCH:
                            gappedSeq.push_back(query.at(qPos++));
                            ++tPos;
                            break;

                        case CigarOperationType::DELETION:
                            gappedSeq.push_back('-');
                            ++tPos;
                            break;

                        case CigarOperationType::INSERTION:  // fall-through
                        case CigarOperationType::SOFT_CLIP:
                            gappedSeq.push_back(query.at(qPos++));
                            break;

                        case CigarOperationType::PADDING:         // fall-through
                        case CigarOperationType::REFERENCE_SKIP:  // .
                        case CigarOperationType::HARD_CLIP:       // .
                        case CigarOperationType::UNKNOWN_OP:
                            throw std::runtime_error("encountered unsupported CIGAR operation");
                    }
                }
            }

            // add gaps after alignment end
            while (tPos < targetSequence.size()) {
                gappedSeq.push_back('-');
                ++tPos;
            }

            // store gapped sequence & insertions for alignment
            alignmentMatrix.push_back(gappedSeq);
            perAlignmentInsertions.push_back(insertions);
        }

        // apply gap-padding across all sequences at insertion sites
        size_t j = 0;
        for (auto& seq : alignmentMatrix) {

            const auto& insertions = perAlignmentInsertions.at(j);
            size_t padsSeen = 0;

            auto insertIter = maxInsertionsLookup.cbegin();
            const auto insertEnd = maxInsertionsLookup.cend();
            for (; insertIter != insertEnd; ++insertIter) {
                auto pos = insertIter->first;
                const auto maxPadsToInsert = insertIter->second;

                // decrement if alignment already has insertions at this position
                auto padsToInsert = maxPadsToInsert;
                const auto existingInserts = insertions.find(pos);
                if (existingInserts != insertions.cend()) padsToInsert -= existingInserts->second;

                // shift position to account for pads already inserted
                pos += padsSeen;

                // insert pads & update padding counter
                seq.insert(pos, std::string(padsToInsert, '-'));
                padsSeen += maxPadsToInsert;
            }
            ++j;
        }

        // sanity check
        std::size_t seqLength = std::string::npos;
        for (const auto& seq : alignmentMatrix) {
            if (seqLength == std::string::npos)
                seqLength = seq.size();
            else
                assert(seqLength == seq.size());
        }

        return alignmentMatrix;
    }

    /// \brief Scan an MSA of sequences for all possible chimeric break-points
    ///        that could explain the query as a composite of the parents
    ///
    /// \param A vector of aligned sequences
    /// \param A string for the name of the Query
    /// \param A string for the name of the first parent
    /// \param A string for the name of the second parent
    /// \param An unsigned int for the row-index of the first parent
    /// \param An unsigned int for the row-index of the second parent
    ///
    /// \return A ChimeraLabel of the highest-scoring possible chimera
    ///
    ChimeraLabel ScorePossibleChimera(const std::vector<std::string>& alignments,
                                      const std::string& queryId, const std::string& parentAId,
                                      const std::string& parentBId, const size_t firstIdx,
                                      const size_t secondIdx) const
    {
        // Extract references to the rows we need to inspect
        const auto& queryRow = alignments.at(0);
        const auto& parentA = alignments.at(firstIdx);
        const auto& parentB = alignments.at(secondIdx);

        // Initialize const-iterators for each row
        auto itQ = queryRow.cbegin();
        auto itA = parentA.cbegin();
        auto itB = parentB.cbegin();

        // Initialize end-point for the const-iterators (since we're traversing
        //    an MSA all 3 iterators will terminate at the same point)
        const auto qEnd = queryRow.cend();

        // Count variables
        size_t rightA = 0;
        size_t rightB = 0;
        size_t rightAbs = 0;
        size_t leftA = 0;
        size_t leftB = 0;
        size_t leftAbs = 0;

        // returns true if input character is a gap
        auto isGap = [](const char c) { return c == '-'; };

        // First we iterate once, counting up the total number of A/B/Abs votes
        //    to initialize the counts for the right-hand segment
        for (; itQ != qEnd; ++itQ, ++itA, ++itB) {

            const char q = *itQ;
            const char a = *itA;
            const char b = *itB;

            // We can't trust our consensus around gaps for low-coverage reads
            //     so we skip them
            if (isGap(q) || isGap(a) || isGap(b)) continue;

            // Identical bases are useless for Chimera Detection - skip
            if (q == a && q == b) continue;

            // If we made it here, count the difference as a vote
            if (a == b) {
                ++rightAbs;
            } else if (q == a) {
                ++rightA;
            } else if (q == b) {
                ++rightB;
            }
        }

        // Re-initialize iterators for each row
        itQ = queryRow.cbegin();
        itA = parentA.cbegin();
        itB = parentB.cbegin();

        // Initialize variables for the maximum chimera
        double maxChimeraScore = 0.0;
        size_t maxChimeraCrossover = 0;
        orientation_t maxChimeraOrientation = NA;

        // And per-iteration variables
        double abScore = 0.0;
        double baScore = 0.0;
        double chimeraScore = 0.0;
        size_t chimeraCrossover = 0;
        orientation_t chimeraOrientation = NA;

        // Second time we iterate, we subtract from the right
        for (; itQ != qEnd; ++itQ, ++itA, ++itB) {

            const char q = *itQ;
            const char a = *itA;
            const char b = *itB;

            // If the Query sequence is not at a gap, increment the
            //     tracking variable for the crossover postion
            if (isGap(q))
                continue;
            else
                ++chimeraCrossover;

            // We can't trust our consensus around gaps for low-coverage reads
            //     so we skip them
            if (isGap(a) || isGap(b)) continue;

            // Identical bases are useless for Chimera Detection - skip
            if (q == a && q == b) continue;

            // If we made it here, count the differences as a vote
            if (a == b) {
                ++leftAbs;
                --rightAbs;
            } else if (q == a) {
                ++leftA;
                --rightA;
            } else if (q == b) {
                ++leftB;
                --rightB;
            }

            // If we've exhausted the right-side votes, there are no more
            //     possible break points
            if (rightA == 0 && rightB == 0) break;

            // If haven't seen any left-side votes, we haven't reached any
            //     possible break points yet
            if (leftA == 0 && leftB == 0) continue;

            // If the Left leans A and the right leans B, test "AAABBB"
            if (leftA > leftB && rightA < rightB) {
                chimeraScore = ScoreBreakPoint(leftA, leftB, leftAbs, rightB, rightA, rightAbs);
                chimeraOrientation = A_B;
            }  // If the Left leans B and the right leans A, test "BBBAAA"
            else if (leftA < leftB && rightA > rightB) {
                chimeraScore = ScoreBreakPoint(leftB, leftA, leftAbs, rightA, rightB, rightAbs);
                chimeraOrientation = B_A;
            }  // If either left or right isn't clear, test both options
            else {
                abScore = ScoreBreakPoint(leftA, leftB, leftAbs, rightB, rightA, rightAbs);
                baScore = ScoreBreakPoint(leftB, leftA, leftAbs, rightA, rightB, rightAbs);
                if (abScore > baScore) {
                    chimeraScore = abScore;
                    chimeraOrientation = A_B;
                } else {
                    chimeraScore = baScore;
                    chimeraOrientation = B_A;
                }
            }

            // Keep the best chimera we've seen so far
            if (chimeraScore > maxChimeraScore) {
                maxChimeraScore = chimeraScore;
                maxChimeraCrossover = chimeraCrossover;
                maxChimeraOrientation = chimeraOrientation;
            }
        }

        // Build and return an appropriately oriented ChimeraLabel
        if (maxChimeraOrientation == A_B)
            return ChimeraLabel(queryId, parentAId, parentBId, maxChimeraCrossover,
                                maxChimeraScore);
        else
            return ChimeraLabel(queryId, parentBId, parentAId, maxChimeraCrossover,
                                maxChimeraScore);
    }

public:  // Public Static Methods
    /// \brief Parse the number of reads supporting a sequence from its Id
    ///
    /// \param The Id to be parsed as a string
    ///
    /// \return The support for that sequence, as an integer
    ///
    static size_t ParseNumReads(const std::string id)
    {
        const auto parts = PacBio::Utility::Split(id, '_');
        const auto numReadsString = parts[3].substr(8);
        const auto numReads = std::stoul(numReadsString);
        return numReads;
    }

    /// \brief Parse the number of reads supporting a vector of sequences from
    /// their Ids
    ///
    /// \param The vector of Ids to be parsed as strings
    ///
    /// \return The support for those sequences, as integers
    ///
    static std::vector<size_t> ParseNumReads(const std::vector<std::string>& ids)
    {
        std::vector<size_t> retval;
        retval.reserve(ids.size());
        for (const auto& id : ids)
            retval.emplace_back(ParseNumReads(id));
        return retval;
    }

private:
    /// \brief Calculates the H-score for a chimeric alignment as per Edgar(2011)
    ///
    /// \param An unsigned int of votes for left-side similarity
    /// \param An unsigned int of votes against left-side similarity
    /// \param An unsigned int of votes neither for or against the left side
    /// \param An unsigned int of votes for right-side similarity
    /// \param An unsigned int of votes against right-side similarity
    /// \param An unsigned int of votes neither for or against the right side
    ///
    /// \return A double representing the strength of the similarity between the
    ///         query sequence and a hypothetical chimera composed of left and
    ///         right segments taken from existing sequences
    ///
    double ScoreBreakPoint(const size_t leftYesVotes, const size_t leftNoVotes,
                           const size_t leftAbsVotes, const size_t rightYesVotes,
                           const size_t rightNoVotes, const size_t rightAbsVotes) const
    {
        // Score the left and right segments independently
        const auto leftScore = ScoreSegment(leftYesVotes, leftNoVotes, leftAbsVotes);
        const auto rightScore = ScoreSegment(rightYesVotes, rightNoVotes, rightAbsVotes);
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
    double ScoreSegment(const size_t yesVotes, const size_t noVotes, const size_t absVotes) const
    {
        static const size_t beta = 4;
        static const double pseudocount = 2.0;
        return yesVotes / (beta * (noVotes + pseudocount) + absVotes);
    }
};

}  // namespace Chimera
}  // namespace PacBio
