// Author: Lance Hepler

#pragma once

#include <memory>
#include <string>
#include <vector>

#include <pbbam/DataSet.h>

namespace PacBio {
namespace IO {

struct SnrInfo
{
    const float A;
    const float C;
    const float G;
    const float T;

    SnrInfo(float a, float c, float g, float t) : A(a), C(c), G(g), T(t) {}
    // comes in TGAC from the instrument, usually
    SnrInfo(float* snr) : SnrInfo(snr[2], snr[3], snr[1], snr[0]) {}
    float Min() const { return std::min(std::min(A, C), std::min(G, T)); }
};

std::vector<std::string> FlattenFofn(const std::vector<std::string>& files);

bool ValidBaseFeatures(const PacBio::BAM::DataSet& ds);

}  // namespace IO
}  // namespace PacBio
