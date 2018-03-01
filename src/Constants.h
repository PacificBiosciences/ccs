// Author: David Seifert

#pragma once

#include <limits>

namespace PacBio {
namespace Consensus {

namespace {

static constexpr const size_t EXTEND_BUFFER_COLUMNS = 8;
static constexpr const int NEG_INT_INF = -std::numeric_limits<int>::infinity();
static constexpr const double NEG_DBL_INF = -std::numeric_limits<double>::infinity();
static constexpr const float NEG_FLOAT_INF = -std::numeric_limits<float>::infinity();

}  // anonymous namespace

}  // namespace Consensus
}  // namespace PacBio
