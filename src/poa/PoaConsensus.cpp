// Author: David Alexander

#include <string>
#include <utility>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>

#include <pacbio/align/AlignConfig.h>
#include <pacbio/denovo/PoaConsensus.h>

namespace PacBio {
namespace Poa {

using namespace PacBio::Align;

AlignConfig DefaultPoaConfig(AlignMode mode)
{
    AlignParams params(3, -5, -4, -4);
    AlignConfig config(params, mode);
    return config;
}

PoaConsensus::PoaConsensus(std::string css, const PoaGraph& g, std::vector<size_t> cssPath)
    : Sequence(std::move(css)), Graph(g), Path(std::move(cssPath))
{
}

PoaConsensus::PoaConsensus(std::string css, const detail::PoaGraphImpl& gi,
                           std::vector<size_t> cssPath)
    : Sequence(std::move(css)), Graph(gi), Path(std::move(cssPath))
{
}

PoaConsensus::~PoaConsensus() = default;
const PoaConsensus* PoaConsensus::FindConsensus(const std::vector<std::string>& reads)
{
    return FindConsensus(reads, DefaultPoaConfig(AlignMode::GLOBAL), -INT_MAX);
}

const PoaConsensus* PoaConsensus::FindConsensus(const std::vector<std::string>& reads,
                                                const AlignConfig& config, int minCoverage)
{
    PoaGraph pg;
    for (const std::string& read : reads) {
        if (read.length() == 0) {
            throw std::invalid_argument("input sequences must have nonzero length.");
        }
        pg.AddRead(read, config);
    }
    return pg.FindConsensus(config, minCoverage);
}

const PoaConsensus* PoaConsensus::FindConsensus(const std::vector<std::string>& reads,
                                                AlignMode mode, int minCoverage)
{
    return FindConsensus(reads, DefaultPoaConfig(mode), minCoverage);
}

std::string PoaConsensus::ToGraphViz(int flags) const { return Graph.ToGraphViz(flags, this); }

void PoaConsensus::WriteGraphVizFile(std::string filename, int flags) const
{
    Graph.WriteGraphVizFile(filename, flags, this);
}

}  // namespace Poa
}  // namespace PacBio
