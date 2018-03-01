// Author: Armin Töpfer

#pragma once

namespace PacBio {
namespace Data {

enum struct StrandType : uint8_t
{
    FORWARD,
    REVERSE,
    UNMAPPED
};
}
}  //::PacBio::Data