#pragma once

#include <ConsensusCore/Poa/PoaGraph.hpp>
#include <ConsensusCore/Interval.hpp>

#include <cstddef>
#include <map>
#include <string>
#include <vector>
#include <utility>

namespace ConsensusCore {
namespace detail {

    // an Anchor represents a point (cssPos, readPos)
    typedef std::pair<size_t, size_t>    SdpAnchor;
    typedef std::vector<SdpAnchor> SdpAnchorVector;

    class PoaGraphImpl;

    using ConsensusCore::PoaGraph;
    using ConsensusCore::Interval;

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
        std::map<PoaGraph::Vertex, Interval> alignableReadIntervalByVertex_;

    public:
        virtual ~SdpRangeFinder();

        void InitRangeFinder(const PoaGraphImpl& poaGraph,
                             const std::vector<PoaGraph::Vertex>& consensusPath,
                             const std::string& consensusSequence,
                             const std::string& readSequence);

        Interval FindAlignableRange(PoaGraph::Vertex v);

    protected:
        virtual SdpAnchorVector FindAnchors(const std::string& consensusSequence,
                                            const std::string& readSequence) const = 0;
    };

}}
