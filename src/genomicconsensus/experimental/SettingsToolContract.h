// Author: Derek Barnett

#pragma once

#include <pbcopper/cli/toolcontract/Config.h>

// clang-format off

namespace PacBio {
namespace GenomicConsensus {
namespace experimental {

inline PacBio::CLI::ToolContract::Config ToolContractConfig()
{
    using Task = PacBio::CLI::ToolContract::Task;

    const std::string id = "genomic_consensus.tasks.gcpp";

    Task tcTask(id);
    tcTask.NumProcessors(Task::MAX_NPROC);
//    tcTask.AddOption();
    tcTask.InputFileTypes({});
    tcTask.OutputFileTypes({});
    return { tcTask };
}

}  // namespace experimental
}  // namespace GenomicConsensus
}  // namespace PacBio

// clang-format on
