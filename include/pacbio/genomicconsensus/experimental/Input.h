// Author: Derek Barnett

#pragma once

#include <string>
#include <vector>

#include <pbbam/BamRecord.h>
#include <pbbam/IndexedFastaReader.h>

#include <pacbio/genomicconsensus/experimental/ReferenceWindow.h>
#include <pacbio/genomicconsensus/experimental/Settings.h>

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

///
/// \brief The Input class
///
class Input
{
public:
    explicit Input(const Settings& settings);

    Input() = delete;
    Input(const Input&) = delete;
    Input(Input&&) = default;
    Input& operator=(const Input&) = delete;
    Input& operator=(Input&&) = default;
    ~Input() = default;

public:
    ///
    /// \brief ReadsInWindow
    /// \param window
    /// \return
    ///
    std::vector<PacBio::BAM::BamRecord> ReadsInWindow(const ReferenceWindow& window) const;

    ///
    /// \brief ReferenceInWindow
    /// \param window
    /// \return
    ///
    std::string ReferenceInWindow(const ReferenceWindow& window) const;

    ///
    /// \brief ReferenceNames
    /// \return
    ///
    std::vector<std::string> ReferenceNames() const;

    ///
    /// \brief ReferenceWindows
    /// \return FASTA references as windows
    ///
    std::vector<ReferenceWindow> ReferenceWindows(bool splitWindows = true) const;

    ///
    /// \brief SequenceLength
    /// \param refName
    /// \return
    ///
    size_t SequenceLength(const std::string& refName) const;

private:
    Settings settings_;
    PacBio::BAM::IndexedFastaReader fasta_;
};

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio
