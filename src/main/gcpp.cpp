// Authors: David Alexander, Lance Hepler, Derek Barnett

#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include <pbcopper/cli/CLI.h>

#include <pacbio/genomicconsensus/experimental/Settings.h>
#include <pacbio/genomicconsensus/experimental/Workflow.h>

int main(int argc, char* argv[])
{
    using namespace PacBio::GenomicConsensus::experimental;

    try {
        return PacBio::CLI::Run(argc, argv, Settings::CreateInterface(), &Workflow::Runner);
    } catch (const std::runtime_error& e) {
        std::cerr << "ERROR: " << e.what();
        exit(EXIT_FAILURE);
    }
}
