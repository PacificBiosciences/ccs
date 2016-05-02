
#include <algorithm>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <boost/foreach.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/optional.hpp>
#include <boost/range/adaptor/reversed.hpp>

#include <pacbio/consensus/poa/RangeFinder.h>

#include "PoaGraphImpl.h"

#define WIDTH 30
#define DEBUG_RANGE_FINDER 0

#if DEBUG_RANGE_FINDER
#include <iostream>
using std::cout;
using std::endl;
#endif  // DEBUG_RANGE_FINDER

namespace PacBio {
namespace Consensus {
namespace detail {

using std::min;
using std::max;
using boost::optional;

static inline bool compareAnchorsOnCssPos(const SdpAnchor& a1, const SdpAnchor& a2)
{
    return a1.first < a2.first;
}

static const SdpAnchor* binarySearchAnchors(const SdpAnchorVector& anchors, size_t cssPosition)
{
    typedef SdpAnchorVector::const_iterator iter_t;
    iter_t found = std::lower_bound(anchors.begin(), anchors.end(), std::make_pair(cssPosition, -1),
                                    compareAnchorsOnCssPos);
    if (found != anchors.end() && (*found).first == cssPosition) {
        return &(*found);
    } else {
        return NULL;
    }
}

typedef std::pair<int, int> Interval;

std::string formatIntervalEndpoint(int i)
{
    if (i == INT_MAX / 2) {
        return "inf";
    } else if (i == -INT_MAX / 2) {
        return "-inf";
    } else {
        return std::to_string(i);
    }
}

std::string formatInterval(const Interval& ival)
{
    return (std::string("[") + formatIntervalEndpoint(ival.first) + std::string(", ") +
            formatIntervalEndpoint(ival.second) + std::string(")"));
}

// Canonical empty interval....
const Interval emptyInterval = Interval(INT_MAX / 2, -INT_MAX / 2);

inline Interval RangeUnion(const Interval& range1, const Interval& range2)
{
    return Interval(min(range1.first, range2.first), max(range1.second, range2.second));
}

inline Interval RangeUnion(const std::vector<Interval>& ranges)
{
    Interval result = emptyInterval;
    for (const Interval& r : ranges) {
        result = RangeUnion(result, r);
    }
    return result;
}

inline Interval next(const Interval& v, int upperBound)
{
    if (v == emptyInterval)
        return emptyInterval;
    else
        return Interval(min(v.first + 1, upperBound), min(v.second + 1, upperBound));
}

inline Interval prev(const Interval& v, int lowerBound = 0)
{
    if (v == emptyInterval)
        return emptyInterval;
    else
        return Interval(max(v.first - 1, lowerBound), max(v.second - 1, lowerBound));
}

SdpRangeFinder::~SdpRangeFinder() {}
void SdpRangeFinder::InitRangeFinder(const PoaGraphImpl& poaGraph,
                                     const std::vector<Vertex>& consensusPath,
                                     const std::string& consensusSequence,
                                     const std::string& readSequence)
{
#if DEBUG_RANGE_FINDER
    poaGraph.WriteGraphVizFile("debug-graph.dot", PoaGraph::VERBOSE_NODES, NULL);
    std::map<Vertex, const SdpAnchor*> anchorByVertex;
#endif
    // Clear prexisting state first!
    alignableReadIntervalByVertex_.clear();

    const int readLength = readSequence.size();
    SdpAnchorVector anchors = FindAnchors(consensusSequence, readSequence);
#if DEBUG_RANGE_FINDER
    std::cout << "Consensus: " << consensusSequence << std::endl;
    std::cout << "Read     : " << readSequence << std::endl;
    std::cout << "RawAnchors length: " << anchors.size() << std::endl;
#endif

    std::map<VD, optional<Interval>> directRanges;
    std::map<VD, Interval> fwdMarks, revMarks;

    std::vector<VD> sortedVertices(num_vertices(poaGraph.g_));
    topological_sort(poaGraph.g_, sortedVertices.rbegin());
    for (const VD v : sortedVertices) {
        directRanges[v] = boost::none;
    }

    // Find the "direct ranges" implied by the anchors between the
    // css and this read.  Possibly null.
    for (size_t cssPos = 0; cssPos < consensusPath.size(); cssPos++) {
        Vertex vExt = consensusPath[cssPos];
        VD v = poaGraph.internalize(vExt);
        const SdpAnchor* anchor = binarySearchAnchors(anchors, cssPos);
        if (anchor != NULL) {
#if DEBUG_RANGE_FINDER
            anchorByVertex[vExt] = anchor;
#endif
            directRanges[v] = Interval(max(int(anchor->second) - WIDTH, 0),
                                       min(int(anchor->second) + WIDTH, readLength));
        } else {
            directRanges[v] = boost::none;
        }
    }

    // Use the direct ranges as a seed and perform a forward recursion,
    // letting a node with null direct range have a range that is the
    // union of the "forward stepped" ranges of its predecessors
    for (const VD v : sortedVertices) {
        Vertex vExt = poaGraph.externalize(v);  // DEBUGGING

        optional<Interval> directRange = directRanges.at(v);
        if (directRange) {
            fwdMarks[v] = directRange.get();
        } else {
            std::vector<Interval> predRangesStepped;
            for (const ED& e : inEdges(v, poaGraph.g_)) {
                VD pred = source(e, poaGraph.g_);
                Vertex predExt = poaGraph.externalize(pred);  // DEBUGGING
                Interval predRangeStepped = next(fwdMarks.at(pred), readLength);
                predRangesStepped.push_back(predRangeStepped);
            }
            Interval fwdInterval = RangeUnion(predRangesStepped);
            fwdMarks[v] = fwdInterval;
        }
    }

    // Do the same thing, but as a backwards recursion
    for (const VD v : boost::adaptors::reverse(sortedVertices)) {
        Vertex vExt = poaGraph.externalize(v);  // DEBUGGING

        optional<Interval> directRange = directRanges.at(v);
        if (directRange) {
            revMarks[v] = directRange.get();
        } else {
            std::vector<Interval> succRangesStepped;
            BOOST_FOREACH (const ED& e, out_edges(v, poaGraph.g_)) {
                VD succ = target(e, poaGraph.g_);
                Vertex succExt = poaGraph.externalize(succ);  // DEBUGGING
                Interval succRangeStepped = prev(revMarks.at(succ), 0);
                succRangesStepped.push_back(succRangeStepped);
            }
            Interval revInterval = RangeUnion(succRangesStepped);
            revMarks[v] = revInterval;
        }
    }

    // take hulls of extents from forward and reverse recursions
    for (const VD v : sortedVertices) {
        Vertex vExt = poaGraph.externalize(v);
        alignableReadIntervalByVertex_[vExt] = RangeUnion(fwdMarks.at(v), revMarks.at(v));
#if DEBUG_RANGE_FINDER
        cout << vExt << "\t";
        if (anchorByVertex.find(vExt) != anchorByVertex.end()) {
            cout << " @  " << anchorByVertex.at(vExt)->second << "\t";
        } else {
            cout << "\t";
        }
        cout << "Fwd mark: " << formatInterval(fwdMarks[v]) << "\t"
             << "Rev mark: " << formatInterval(revMarks[v]) << "\t"
             << "Range: " << formatInterval(alignableReadIntervalByVertex_[vExt]) << endl;
#endif
    }
}

Interval SdpRangeFinder::FindAlignableRange(Vertex v)
{
    return alignableReadIntervalByVertex_.at(v);
}

}  // namespace detail
}  // namespace Consensus
}  // namespace PacBio
