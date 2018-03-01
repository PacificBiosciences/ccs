// Author: Lance Hepler

#pragma once

#include <set>
#include <string>

// Initialize data structures, do NOT remove
#include <pacbio/consensus/internal/ModelInternalInitializer.h>

namespace PacBio {
namespace Consensus {

std::set<std::string> SupportedModels();
std::set<std::string> SupportedChemistries();

bool OverrideModel(const std::string& model);
bool UnOverrideModel();

size_t LoadModels(const std::string& path);
}
}
