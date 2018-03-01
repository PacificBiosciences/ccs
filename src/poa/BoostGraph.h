// Author: David Alexander

#pragma once

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>

namespace boost {
enum vertex_info_t
{
    vertex_info = 424
};  // a unique #
BOOST_INSTALL_PROPERTY(vertex, info);
}  // namespace boost

namespace PacBio {
namespace Poa {
namespace detail {

/* using boost::adjacency_list; */
/* using boost::graph_traits; */
/* using boost::property; */
/* using boost::property_map; */

using namespace boost;

struct PoaNode;

// BGL is intimidating, and it *deserves* your hatred.  But it's
// the only game in town!
using vertex_property_t = property<vertex_info_t, PoaNode, property<vertex_index_t, size_t> >;
using BoostGraph = adjacency_list<setS, listS, bidirectionalS, vertex_property_t>;

// Descriptor types used internally
using ED = graph_traits<BoostGraph>::edge_descriptor;
using VD = graph_traits<BoostGraph>::vertex_descriptor;

using VertexInfoMap = property_map<BoostGraph, vertex_info_t>::type;
using index_map_t = property_map<BoostGraph, vertex_index_t>::type;
static const VD null_vertex = graph_traits<BoostGraph>::null_vertex();
}
}
}
