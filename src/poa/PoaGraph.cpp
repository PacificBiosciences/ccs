// Author: David Alexander
// (Based on the original "Partial Order Aligner" by Lee, Grasso, and Sharlow,
//  and an implementation in C# by Patrick Marks)

#include <pacbio/denovo/PoaGraph.h>

#include "PoaGraphImpl.h"

using namespace PacBio::Align;

namespace PacBio {
namespace Poa {

// forward declaration
struct PoaConsensus;

//
// PIMPL idiom delegation
//
void PoaGraph::AddRead(const std::string& sequence, const AlignConfig& config,
                       detail::SdpRangeFinder* rangeFinder, std::vector<Vertex>* readPathOutput)
{
    impl->AddRead(sequence, config, rangeFinder, readPathOutput);
}

void PoaGraph::AddFirstRead(const std::string& sequence, std::vector<Vertex>* readPathOutput)
{
    impl->AddFirstRead(sequence, readPathOutput);
}

PoaAlignmentMatrix* PoaGraph::TryAddRead(const std::string& sequence, const AlignConfig& config,
                                         detail::SdpRangeFinder* rangeFinder) const
{
    return impl->TryAddRead(sequence, config, rangeFinder);
}

void PoaGraph::CommitAdd(PoaAlignmentMatrix* mat, std::vector<Vertex>* readPathOutput)
{
    impl->CommitAdd(mat, readPathOutput);
}

void PoaGraph::PruneGraph(const int minCoverage) { impl->PruneGraph(minCoverage); }

size_t PoaGraph::NumReads() const { return impl->NumReads(); }

const PoaConsensus* PoaGraph::FindConsensus(const AlignConfig& config, int minCoverage) const
{
    return impl->FindConsensus(config, minCoverage);
}

string PoaGraph::ToGraphViz(int flags, const PoaConsensus* pc) const
{
    return impl->ToGraphViz(flags, pc);
}

void PoaGraph::WriteGraphVizFile(const string& filename, int flags, const PoaConsensus* pc) const
{
    impl->WriteGraphVizFile(filename, flags, pc);
}

void PoaGraph::WriteGraphCsvFile(const string& filename) const
{
    impl->WriteGraphCsvFile(filename);
}

PoaGraph::PoaGraph() { impl = new detail::PoaGraphImpl(); }

PoaGraph::PoaGraph(const PoaGraph& other) { impl = new detail::PoaGraphImpl(*other.impl); }

PoaGraph::PoaGraph(const detail::PoaGraphImpl& o) { impl = new detail::PoaGraphImpl(o); }

PoaGraph::~PoaGraph() { delete impl; }

}  // namespace Poa
}  // namespace PacBio
