// Authors: Lance Hepler, Brett Bowman

#pragma once

#include <vector>

#include <seqan/basic.h>

#include <pacbio/data/internal/BaseEncoding.h>

namespace PacBio {
namespace Align {

template <typename TShape>
class HomopolymerHasher
{

public:  // Structors
    // Default constructor
    HomopolymerHasher(TShape& shape) { Initialize(shape); }

    // Move constructor
    HomopolymerHasher(HomopolymerHasher&& src) = delete;
    // Copy constructor
    HomopolymerHasher(const HomopolymerHasher& src) = delete;
    // Move assignment operator
    HomopolymerHasher& operator=(HomopolymerHasher&& rhs) = delete;
    // Copy assignment operator
    HomopolymerHasher& operator=(const HomopolymerHasher& rhs) = delete;
    // Destructor
    ~HomopolymerHasher() = default;

private:  // Internal class methods
    /// Initialize the hasher object by creating all possible
    /// homopolymer templates that fill the shape used by the
    /// index and saves them for later used by
    ///
    /// \param  shape  The Q-gram shape used in the index
    void Initialize(const TShape& shape)
    {
        using namespace seqan;

        for (size_t i = 0; i < 4; i++) {
            DnaString s = std::string(length(shape), Data::detail::NCBI2naToASCIIImpl(i));
            hashes[i] = seqan::hash(shape, begin(s));
        }
    }

public:  // Main Function
    /// Given the hash of a Q-Gram, check whether it matches the
    /// pattern of one of the stored hashes for known homopolymers.
    /// If so return True, otherwise False.
    ///
    /// \param  h  The Q-Gram hash of the query sequence
    inline bool operator()(const unsigned h) const
    {
        if (h == hashes[0] || h == hashes[1] || h == hashes[2] || h == hashes[3]) return true;

        return false;
    }

private:  // Class variables
    unsigned hashes[4];

};  // class HomopolymerHasher

}  // namespace Align
}  // namespace PacBio
