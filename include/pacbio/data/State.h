// Author: Armin Töpfer

#pragma once

#include <iostream>

namespace PacBio {
namespace Data {

enum struct State : uint8_t
{
    VALID = 0,
    ALPHA_BETA_MISMATCH,
    POOR_ZSCORE,
    TEMPLATE_TOO_SMALL,
    MANUALLY_RELEASED,
    ILLEGAL_BASE,
    ILLEGAL_PW,
    INVALID,

    SIZE
};

static constexpr const char* StateName[] = {"VALID",
                                            "ALPHA/BETA MISMATCH",
                                            "POOR Z-SCORE",
                                            "TEMPLATE TOO SMALL",
                                            "MANUALLY RELEASED",
                                            "ILLEGAL BASE",
                                            "ILLEGAL PULSEWIDTH",
                                            "INVALID"};

inline std::ostream& operator<<(std::ostream& os, State result)
{
    os << PacBio::Data::StateName[static_cast<size_t>(result)];
    return os;
}
}
}  //::PacBio::Data