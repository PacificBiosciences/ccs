// Author: David Alexander

#pragma once

#include <climits>
#include <string>
#include <utility>
#include <vector>

#include <pacbio/align/AlignConfig.h>

namespace PacBio {
namespace Poa {

// fwd decls
namespace detail {

class PoaGraphImpl;
class SdpRangeFinder;

}  // namespace detail

struct PoaConsensus;

class PoaAlignmentMatrix
{
public:
    virtual ~PoaAlignmentMatrix(){};
    virtual float Score() const = 0;
    virtual size_t NumRows() const = 0;
    virtual size_t NumCols() const = 0;
    virtual void Print() const = 0;
};

/// \brief An object representing a Poa (partial-order alignment) graph
class PoaGraph
{
public:
    typedef size_t Vertex;
    typedef size_t ReadId;

    static const Vertex NullVertex = (Vertex)-1;

public:  // Flags enums for specifying GraphViz output features
    enum
    {
        COLOR_NODES = 0x1,
        VERBOSE_NODES = 0x2
    };

public:
    PoaGraph();
    PoaGraph(const PoaGraph& other);
    PoaGraph(const detail::PoaGraphImpl& o);  // NB: this performs a copy
    ~PoaGraph();

    //
    // Easy API
    //
    void AddRead(const std::string& sequence, const PacBio::Align::AlignConfig& config,
                 detail::SdpRangeFinder* rangeFinder = NULL,
                 std::vector<Vertex>* readPathOutput = NULL);

    //
    // API for more control
    //
    void AddFirstRead(const std::string& sequence, std::vector<Vertex>* readPathOutput = NULL);

    PoaAlignmentMatrix* TryAddRead(const std::string& sequence,
                                   const PacBio::Align::AlignConfig& config,
                                   detail::SdpRangeFinder* rangeFinder = NULL) const;

    void CommitAdd(PoaAlignmentMatrix* mat, std::vector<Vertex>* readPathOutput = NULL);

    void PruneGraph(const int minCoverage);

    // ----------

    size_t NumReads() const;

    std::string ToGraphViz(int flags = 0, const PoaConsensus* pc = NULL) const;

    void WriteGraphVizFile(const std::string& filename, int flags = 0,
                           const PoaConsensus* pc = NULL) const;

    void WriteGraphCsvFile(const std::string& filename) const;

    const PoaConsensus* FindConsensus(const PacBio::Align::AlignConfig& config,
                                      int minCoverage = -INT_MAX) const;

private:
    detail::PoaGraphImpl* impl;
};

}  // namespace Poa
}  // namespace PacBio
