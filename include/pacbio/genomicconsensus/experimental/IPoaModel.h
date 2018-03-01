// Author: Derek Barnett

#pragma once

#include <string>
#include <vector>

#include <pbbam/BamRecord.h>

#include <pacbio/genomicconsensus/experimental/IConsensusModel.h>
#include <pacbio/genomicconsensus/experimental/Variant.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

struct Consensus;
class Input;
struct ReferenceWindow;
struct Settings;
struct WorkChunk;

///
/// \brief The IPoaModel class
///
class IPoaModel : public IConsensusModel
{
public:
    virtual ~IPoaModel();

    /// \brief ProcessChunk
    ///
    ///
    /// \param chunk
    /// \param settings
    /// \return
    ///
    WindowResult ProcessChunk(const WorkChunk& chunk, const Settings& settings);

    ///
    /// \brief ConsensusAndVariantsFromWindow
    /// \param input
    /// \param reads
    /// \param window
    /// \param refSeq
    /// \param settings
    /// \return
    ///
    virtual WindowResult ConsensusAndVariantsFromWindow(
        const Input& input, const std::vector<PacBio::BAM::BamRecord>& reads,
        const ReferenceWindow& window, const std::string& refSeq,
        const Settings& settings) const = 0;

    // ----------------------------------------------------------- //
    // TODO (DB): stateless methods, move out for easier testing   //
    //            ?? IPoaUtils.h ??                               //

    void AnnotateVariants(std::vector<Variant>* const variants,
                          const std::vector<PacBio::BAM::BamRecord>& reads) const;

    void ClipReadsToWindow(std::vector<PacBio::BAM::BamRecord>* const reads,
                           const ReferenceWindow& subWindow) const;

    ReferenceWindow EnlargedWindow(const ReferenceWindow& window, const size_t maxSeqLength,
                                   const size_t overlap) const;

    Consensus RestrictedConsensus(const Consensus& enlargedCss, const std::string& refSeq,
                                  const ReferenceWindow& originalWindow) const;

    std::vector<Variant> RestrictedVariants(const std::vector<Variant>& enlargedVariants,
                                            const ReferenceWindow& originalWindow) const;

    WindowResult ResultForWindow(const ReferenceWindow& window, const std::string& refSeq,
                                 const Settings& settings) const;

    //                                                             //
    // ----------------------------------------------------------- //

protected:
    IPoaModel() = default;
};

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
