// Author: David Seifert

#pragma once

#include <array>
#include <cassert>
#include <cstdint>
#include <stdexcept>

#include <pacbio/UnanimityConfig.h>

#include "ConversionFunctions.h"

namespace PacBio {
namespace Data {
namespace detail {

// Unanimity uses 3 different formats internally for representing
// bases: ASCII, NCBI2na and NCBI4na
//
//   1. ASCII is the standard format and encoding known
//      from FASTA files, that is the base A is encoded as 'A'
//      with ordinal value 65. While ASCII values are nice for
//      strings and human interaction, they are cumbersome for
//      algorithms.
//
//   2. The NCBI2na format is a 2-bit format that does not
//      allow for encoding ambiguous bases:
//
//        Bit2 Bit1 | Dec | State
//        -----------------------
//         0    0   |  0  | A
//         0    1   |  1  | C
//         1    0   |  2  | G
//         1    1   |  3  | T
//
//   3. The NCBI4na format is a 4-bit nibble format that can
//      encode ambiguous bases and is used extensively for
//      diploid usecases:
//
//         T   G   C   A  |  Dec | State
//        -----------------------------------------------
//         0   0   0   0  |  0   | INVALID (gap in NCBI4na)
//         0   0   0   1  |  1   | A
//         0   0   1   0  |  2   | C
//         0   0   1   1  |  3   | M (A/C)
//         0   1   0   0  |  4   | G
//         0   1   0   1  |  5   | R (A/G)
//         0   1   1   0  |  6   | S (C/G)
//         0   1   1   1  |  7   | V (A/C/G)
//         1   0   0   0  |  8   | T
//         1   0   0   1  |  9   | W (A/T)
//         1   0   1   0  |  10  | Y (C/T)
//         1   0   1   1  |  11  | H (A/C/T)
//         1   1   0   0  |  12  | K (G/T)
//         1   1   0   1  |  13  | D (A/G/T)
//         1   1   1   0  |  14  | B (C/G/T)
//         1   1   1   1  |  15  | N (A/C/G/T)
//
// Reference:
//   https://www.ncbi.nlm.nih.gov/IEB/ToolBox/SDKDOCS/BIOSEQ.HTML
//
// IMPORTANT NOTICE
// It should be emphasized here again, that all functions in this file
// are an *implementation detail* of Unanimity. Unanimity does not
// leak ambiguous bases across public API boundaries. You cannot rely
// on any interfaces in this file when consuming Unanimity. You have
// been warned.

// fwd declaration for conversion constructor
class NCBI4na;

class NCBI2na
{
public:
    static inline UNANIMITY_CONSTEXPR NCBI2na FromASCII(const char base) { return NCBI2na{base}; }

    static inline UNANIMITY_CONSTEXPR NCBI2na FromRaw(const uint8_t raw) { return NCBI2na{raw}; }

public:
    ~NCBI2na() = default;

    inline UNANIMITY_CONSTEXPR NCBI2na(const NCBI2na&) = default;
    inline UNANIMITY_CONSTEXPR NCBI2na(NCBI2na&&) = default;

    inline UNANIMITY_CONSTEXPR NCBI2na& operator=(const NCBI2na&) = default;
    inline UNANIMITY_CONSTEXPR NCBI2na& operator=(NCBI2na&&) = default;

    friend class NCBI4na;

public:
    // Accessors
    inline UNANIMITY_CONSTEXPR const uint8_t& Data() const { return data_; }

    inline UNANIMITY_CONSTEXPR NCBI4na GetNCBI4na() const;

    inline UNANIMITY_CONSTEXPR char GetASCII() const { return NCBI2naToASCIIImpl(data_); }

    inline UNANIMITY_CONSTEXPR bool IsValid() const { return (data_ < 4); }

    inline UNANIMITY_CONSTEXPR bool IsEqual(const NCBI2na& rhs) const
    {
        return (data_ == rhs.data_);
    }

private:
    uint8_t data_;

private:
    explicit inline UNANIMITY_CONSTEXPR NCBI2na(const char base) : data_{ASCIIToNCBI2naImpl(base)}
    {
    }
    explicit inline UNANIMITY_CONSTEXPR NCBI2na(const NCBI4na base);
    explicit inline UNANIMITY_CONSTEXPR NCBI2na(const uint8_t raw) : data_{raw} {}
};

class NCBI4na
{
public:
    static inline UNANIMITY_CONSTEXPR NCBI4na FromASCII(const char base,
                                                        const bool checkValid = true)
    {
        return NCBI4na{base, checkValid};
    }

public:
    ~NCBI4na() = default;

    inline UNANIMITY_CONSTEXPR NCBI4na(const NCBI4na&) = default;
    inline UNANIMITY_CONSTEXPR NCBI4na(NCBI4na&&) = default;

    inline UNANIMITY_CONSTEXPR NCBI4na& operator=(const NCBI4na&) = default;
    inline UNANIMITY_CONSTEXPR NCBI4na& operator=(NCBI4na&&) = default;

    friend class NCBI2na;

public:
    // Accessors
    inline UNANIMITY_CONSTEXPR const uint8_t& Data() const { return data_; }

    inline UNANIMITY_CONSTEXPR NCBI2na GetNCBI2na() const
    {
        return NCBI2na{NCBI4naToNCBI2naImpl(data_)};
    }

    inline UNANIMITY_CONSTEXPR char GetASCII() const { return NCBI4naToASCIIImpl(data_); }

    inline UNANIMITY_CONSTEXPR bool Overlap(const NCBI4na& rhs) const
    {
        return (data_ & rhs.data_);
    }

    inline UNANIMITY_CONSTEXPR bool IsValid() const { return ((data_ > 0) && (data_ < 16)); }

    inline UNANIMITY_CONSTEXPR uint8_t NumContainedBases() const { return numSetBitsImpl(data_); }

    inline UNANIMITY_CONSTEXPR bool Contains(const NCBI2na& base) const
    {
        const NCBI4na testBase = base.GetNCBI4na();
        return Overlap(testBase);
    }

    // A/C/G/T are pure bases
    inline UNANIMITY_CONSTEXPR bool IsPure() const { return (NumContainedBases() == 1); }

    inline UNANIMITY_CONSTEXPR bool IsAmbig() const { return (NumContainedBases() > 1); }

    inline UNANIMITY_CONSTEXPR bool IsEqual(const NCBI4na& rhs) const
    {
        return (data_ == rhs.data_);
    }

private:
    uint8_t data_;

private:
    explicit inline UNANIMITY_CONSTEXPR NCBI4na(const char base, const bool checkValid)
        : data_{ASCIIToNCBI4naImpl(base, checkValid)}
    {
    }
    explicit inline UNANIMITY_CONSTEXPR NCBI4na(const NCBI2na base)
        : data_{NCBI2naToNCBI4naImpl(base.data_)}
    {
    }
    explicit inline UNANIMITY_CONSTEXPR NCBI4na(const uint8_t raw) : data_{raw} {}
};

inline UNANIMITY_CONSTEXPR NCBI4na NCBI2na::GetNCBI4na() const
{
    return NCBI4na{NCBI2naToNCBI4naImpl(data_)};
}

inline UNANIMITY_CONSTEXPR NCBI2na::NCBI2na(const NCBI4na base)
    : data_{NCBI4naToNCBI2naImpl(base.data_)}
{
}

}  // namespace detail
}  // namespace Data
}  // namespace PacBio
