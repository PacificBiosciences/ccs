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
namespace Consensus {
namespace detail {

/* using boost::adjacency_list; */
/* using boost::graph_traits; */
/* using boost::property; */
/* using boost::property_map; */

using namespace boost;

struct PoaNode;

// BGL is intimidating, and it *deserves* your hatred.  But it's
// the only game in town!
typedef property<vertex_info_t, PoaNode, property<vertex_index_t, size_t> > vertex_property_t;
typedef adjacency_list<setS, listS, bidirectionalS, vertex_property_t> BoostGraph;

// Descriptor types used internally
typedef graph_traits<BoostGraph>::edge_descriptor ED;
typedef graph_traits<BoostGraph>::vertex_descriptor VD;

typedef property_map<BoostGraph, vertex_info_t>::type VertexInfoMap;
typedef property_map<BoostGraph, vertex_index_t>::type index_map_t;
static const VD null_vertex = graph_traits<BoostGraph>::null_vertex();
}
}
}
