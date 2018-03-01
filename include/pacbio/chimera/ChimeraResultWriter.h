// Author: Brett Bowman

#pragma once

#include <vector>

#include "ChimeraLabel.h"

namespace PacBio {
namespace Chimera {

///
/// \brief Write out to file a series ChimeraLabels in human-readable form
///
class ChimeraResultWriter
{
public:
    ChimeraResultWriter(const std::string &csvFile) : csvWriter_(csvFile)
    {
        ValidateCsvFile(csvFile);
        csvWriter_ << "SequenceId,IsChimera,ChimeraScore,"
                   << "ParentSequenceA,ParentSequenceB,CrossoverPosition" << std::endl;
    }

    void WriteResult(const ChimeraLabel &label) { csvWriter_ << label << std::endl; }
    void WriteResults(const std::vector<ChimeraLabel> &labels)
    {
        for (const auto &label : labels)
            WriteResult(label);
    }

private:
    static void ValidateCsvFile(const std::string &filename)
    {
        const std::string ext = filename.substr(filename.find_last_of(".") + 1);

        if (ext != "csv")
            throw std::invalid_argument("invalid sequence file type, only CSVs supported");
    }

private:
    std::ofstream csvWriter_;
};

}  // namespace Chimera
}  // namespace PacBio
