// Author: David Alexander

#pragma once

#include <string>
#include <vector>

bool LoadFastaSequences(std::string fastaFname, std::vector<std::string>& ids,
                        std::vector<std::string>& sequences);
