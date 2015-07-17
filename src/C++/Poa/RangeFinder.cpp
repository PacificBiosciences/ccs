
#include <ConsensusCore/Poa/RangeFinder.hpp>
#include <ConsensusCore/Interval.hpp>

#include <algorithm>
#include <boost/graph/topological_sort.hpp>
#include <boost/optional.hpp>
#include <map>
#include <string>
#include <vector>
#include <utility>

#include "PoaGraphImpl.hpp"

#define WIDTH 30
#define DEBUG_RANGE_FINDER 0

#if DEBUG_RANGE_FINDER
#include <iostream>
using std::cout;
using std::endl;
#endif  // DEBUG_RANGE_FINDER


namespace ConsensusCore {
namespace detail {

    using std::min;
    using std::max;
    using std::make_pair;
    using boost::optional;

    static inline bool compareAnchorsOnCssPos(const SdpAnchor& a1, const SdpAnchor& a2)
    {
        return a1.first < a2.first;
    }

    static const SdpAnchor* binarySearchAnchors(const SdpAnchorVector& anchors, size_t cssPosition)
    {
        typedef SdpAnchorVector::const_iterator iter_t;
        iter_t found = std::lower_bound(anchors.begin(), anchors.end(),
                                        make_pair(cssPosition, -1), compareAnchorsOnCssPos);
        if (found != anchors.end() && (*found).first == cssPosition)
        { return &(*found); }
        else
        { return NULL; }
    }


    static inline Interval next(const Interval& v, int upperBound)
    {
        return Interval(min(v.Begin + 1, upperBound),
                        min(v.End + 1,   upperBound));
    }

    static inline Interval prev(const Interval& v, int lowerBound=0)
    {
        return Interval(max(v.Begin - 1, lowerBound),
                        max(v.End - 1,   lowerBound));
    }

    inline Interval
    rangeIntersection(const Interval& range1, const Interval& range2)
    {
        return Interval(max(range1.Begin, range2.Begin),
                        min(range1.End,   range2.End));
    }

    SdpRangeFinder::~SdpRangeFinder() {}

    void
    SdpRangeFinder::InitRangeFinder(const PoaGraphImpl& poaGraph,
                                    const std::vector<Vertex>& consensusPath,
                                    const std::string& consensusSequence,
                                    const std::string& readSequence)
    {
#if DEBUG_RANGE_FINDER
        poaGraph.WriteGraphVizFile("debug-graph.dot", PoaGraph::VERBOSE_NODES, NULL);
#endif
        // Clear prexisting state first!
        alignableReadIntervalByVertex_.clear();

        const int readLength = readSequence.size();

        SdpAnchorVector anchors = FindAnchors(consensusSequence, readSequence);

        std::map<VD, optional<Interval> > directRanges;
        std::map<VD, Interval> fwdMarks, revMarks;

        std::vector<VD> sortedVertices(num_vertices(poaGraph.g_));
        topological_sort(poaGraph.g_, sortedVertices.rbegin());
        foreach (VD v, sortedVertices)
        {
            directRanges[v] = boost::none;
        }

        // Find the "direct ranges" implied by the anchors between the
        // css and this read.  Possibly null.
        for (size_t cssPos = 0; cssPos < consensusPath.size(); cssPos++)
        {
            Vertex vExt = consensusPath[cssPos];
            VD v = poaGraph.internalize(vExt);
            const SdpAnchor* anchor = binarySearchAnchors(anchors, cssPos);
            if (anchor != NULL) {
#if DEBUG_RANGE_FINDER
                cout << "Anchor: " << anchor->first << "-" << anchor->second
                     <<  " (Vertex " << vExt << ")" << endl;
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
        foreach (VD v, sortedVertices)
        {
            optional<Interval> directRange = directRanges.at(v);
            if (directRange) {
                fwdMarks[v] = directRange.get();
            } else {
                std::vector<Interval> predRangesStepped;
                foreach (ED e, inEdges(v, poaGraph.g_))
                {
                    VD pred = source(e, poaGraph.g_);
                    Interval predRangeStepped = next(fwdMarks.at(pred), readLength);
                    predRangesStepped.push_back(predRangeStepped);
                }
                fwdMarks[v] = RangeUnion(predRangesStepped);
            }
        }

        // Do the same thing, but as a backwards recursion
        foreach (VD v, make_pair(sortedVertices.rbegin(), sortedVertices.rend()))
        {
            optional<Interval> directRange = directRanges.at(v);
            if (directRange) {
                revMarks[v] = directRange.get();
            } else {
                std::vector<Interval> succRangesStepped;
                foreach (ED e, out_edges(v, poaGraph.g_))
                {
                    VD succ = target(e, poaGraph.g_);
                    Interval succRangeStepped = prev(revMarks.at(succ), 0);
                    succRangesStepped.push_back(succRangeStepped);
                }
                revMarks[v] = RangeUnion(succRangesStepped);
            }
        }

        // take hulls of extents from forward and reverse recursions
        foreach (VD v, sortedVertices)
        {
            Vertex vExt = poaGraph.externalize(v);
            alignableReadIntervalByVertex_[vExt] = RangeUnion(fwdMarks.at(v), revMarks.at(v));
#if DEBUG_RANGE_FINDER
            cout << vExt << " range = ["
                 << alignableReadIntervalByVertex_[v].Begin  << ", "
                 << alignableReadIntervalByVertex_[v].End << ")" << endl;
#endif
        }
    }

    Interval
    SdpRangeFinder::FindAlignableRange(Vertex v)
    {
        return alignableReadIntervalByVertex_.at(v);
    }

}}
