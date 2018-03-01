// Author: David Alexander

#pragma once

#include <cstdint>

namespace PacBio {
namespace Align {

//
// Scoring params for Needleman-Wunsch or Smith-Waterman style aligners
//
struct AlignParams
{
    int Match;
    int Mismatch;
    int Insert;
    int Delete;

    AlignParams(int match, int mismatch, int insert, int delete_);

    // Edit distance params
    static AlignParams Default();
};

enum struct AlignMode : uint8_t
{
    GLOBAL = 0,      // Global in both sequences
    SEMIGLOBAL = 1,  // Global in query, local in target
    LOCAL = 2        // Local in both sequences
};

struct AlignConfig
{
    AlignParams Params;
    AlignMode Mode;

    AlignConfig(AlignParams params, AlignMode mode);

    // Default corresponds to global alignment mode, edit distance params
    static AlignConfig Default();
};

}  // namespace Align
}  // namespace PacBio
