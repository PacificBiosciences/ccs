// Author: Brett Bowman

#include <seqan/seq_io.h>

#include "pacbio/chimera/ChimeraLabeler.h"
#include "pacbio/chimera/ChimeraResultWriter.h"
#include "pbbam/FastaReader.h"
#include "pbbam/FastqReader.h"

using PacBio::BAM::FastqReader;
using PacBio::BAM::FastaReader;
using namespace PacBio::Chimera;

int main(int argc, char const** argv)
{
    using namespace seqan;

    if (argc != 2) return 1;  // Invalid number of arguments.

    // Get the input file
    std::string inputFile(argv[1]);

    // Parse the records
    std::vector<std::string> ids;
    std::vector<std::string> seqs;
    for (const auto& seq : FastqReader::ReadAll(inputFile)) {
        ids.push_back(seq.Name());
        seqs.push_back(seq.Bases());
    }

    // Label the Records
    ChimeraLabeler chimeraLabeler(1.0f, 100, true);
    auto labels = chimeraLabeler.LabelChimeras(ids, seqs);

    // Display the results
    ChimeraResultWriter csvWriter("temp.csv");
    csvWriter.WriteResults(labels);

    return 0;
}
