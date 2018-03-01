// Author: Armin Töpfer

#pragma once

#include <string>
#include <tuple>

namespace PacBio {

std::string UnanimityChangelog();
std::string UnanimityGitSha1();
std::string UnanimityVersion();
std::tuple<int, int, int> UnanimityVersionTriple();

}  // ::PacBio
