// Author: Lance Hepler

#pragma once

#include <boost/property_tree/ptree.hpp>

namespace PacBio {
namespace Consensus {
template <size_t I>
void ReadMatrix(double (&mat)[I], const boost::property_tree::ptree& pt)
{
    if (pt.size() != I) throw std::invalid_argument("bad size (1D)");
    size_t i = 0;
    for (const auto& item : pt) {
        mat[i] = item.second.get_value<double>();
        ++i;
    }
}

template <size_t I, size_t J>
void ReadMatrix(double (&mat)[I][J], const boost::property_tree::ptree& pt)
{
    if (pt.size() != I) throw std::invalid_argument("bad size (2D)");
    size_t i = 0;
    for (const auto& item : pt) {
        ReadMatrix<J>(mat[i], item.second);
        ++i;
    }
}

template <size_t I, size_t J, size_t K>
void ReadMatrix(double (&mat)[I][J][K], const boost::property_tree::ptree& pt)
{
    if (pt.size() != I) throw std::invalid_argument("bad size (3D)");
    size_t i = 0;
    for (const auto& item : pt) {
        ReadMatrix<J, K>(mat[i], item.second);
        ++i;
    }
}
}
}
