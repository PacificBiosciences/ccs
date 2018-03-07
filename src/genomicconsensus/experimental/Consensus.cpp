// Author: Derek Barnett

#include <pacbio/genomicconsensus/experimental/Consensus.h>

#include <stdexcept>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

Consensus Consensus::Join(std::vector<Consensus> subconsensi)
{
    //
    // [Consensus] -> Consensus
    //
    // String together all the consensus objects into a single consensus.
    // Will throw if the reference windows are not contiguous.
    //

    if (subconsensi.empty()) throw std::runtime_error("cannot join empty Consensus chunk list");
    std::sort(subconsensi.begin(), subconsensi.end());
    std::vector<ReferenceWindow> windows;
    for (const auto& c : subconsensi)
        windows.push_back(c.window);
    if (!AreContiguous(windows)) throw std::runtime_error("Consensus chunks must be contiguous");

    std::string joinedSeq;
    std::vector<uint8_t> joinedConfidence;
    for (const auto& c : subconsensi) {
        joinedSeq += c.sequence;
        joinedConfidence.insert(joinedConfidence.end(), c.confidence.begin(), c.confidence.end());
    }

    return Consensus{
        ReferenceWindow{subconsensi.front().window.name,
                        {subconsensi.front().window.Start(), subconsensi.back().window.End()}},
        std::move(joinedSeq), std::move(joinedConfidence)};
}

Consensus Consensus::NoCallConsensus(const NoCallStyle style, const ReferenceWindow& window,
                                     const std::string& refSeq)
{
    const auto length = refSeq.size();
    switch (style) {
        case (NoCallStyle::NO_CALL): {
            return Consensus{window, std::string(length, 'N'), std::vector<uint8_t>(length, 0)};
        }
        case (NoCallStyle::REFERENCE): {
            return Consensus{window, refSeq, std::vector<uint8_t>(length, 0)};
        }
        case (NoCallStyle::LOWERCASE_REFERENCE): {
            std::string seq;
            seq.reserve(refSeq.size());
            std::transform(refSeq.begin(), refSeq.end(), std::back_inserter(seq),
                           [](const char c) { return std::tolower(c); });
            return Consensus{window, std::move(seq), std::vector<uint8_t>(length, 0)};
        }
        default: {
            // Silence -Wreturn-type
            throw std::runtime_error{"Unknown call style!"};
        }
    }
}

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
