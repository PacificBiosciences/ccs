#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <pacbio/consensus/poa/PoaGraph.h>

namespace PacBio {
namespace Consensus {
namespace detail {

// an Anchor represents a point (cssPos, readPos)
typedef std::pair<size_t, size_t> SdpAnchor;
typedef std::vector<SdpAnchor> SdpAnchorVector;

class PoaGraphImpl;

// using PacBio::Consensus::PoaGraph;

//
// SdpRangeFinder objects are responsible for identifying the range
// of read positions that we should seek to align to a POA vertex;
// this implementation uses SDP to identify fairly narrow bands,
// enabling sparse memory usage.
//
// This is an abstract class that will be inherited in a client
// library that has access to an SDP method.
//
// RangeFinder state goes away on next call to InitRangeFinder.  We could
// have dealt with this using a factory pattern but bleh.
//
class SdpRangeFinder
{
private:
    std::map<PoaGraph::Vertex, std::pair<int, int>> alignableReadIntervalByVertex_;

public:
    virtual ~SdpRangeFinder();

    void InitRangeFinder(const PoaGraphImpl& poaGraph,
                         const std::vector<PoaGraph::Vertex>& consensusPath,
                         const std::string& consensusSequence, const std::string& readSequence);

    // TODO: write contract
    std::pair<int, int> FindAlignableRange(PoaGraph::Vertex v);

protected:
    // TODO: write contract
    virtual SdpAnchorVector FindAnchors(const std::string& consensusSequence,
                                        const std::string& readSequence) const = 0;
};

}  // namespace detail
}  // namespace Consensus
}  // namespace PacBio
