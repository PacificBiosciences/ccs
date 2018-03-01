// Author: Lance Hepler

#include <climits>
#include <cstdlib>
#include <fstream>

#include <sys/stat.h>

#include <boost/algorithm/string.hpp>

#include <pacbio/io/Utility.h>

using std::string;
using std::vector;

using std::ifstream;

namespace PacBio {
namespace IO {
void FlattenFofn(vector<string>& res, const string& file)
{
    using boost::algorithm::iends_with;
    using boost::algorithm::trim;

    if (iends_with(file, ".fofn")) {
        ifstream fofn(file);
        string line;
        while (getline(fofn, line)) {
            trim(line);
            FlattenFofn(res, line);
        }
    } else if (iends_with(file, ".bam"))
        res.push_back(file);
    else
        throw std::invalid_argument("not a .fofn or .bam file!");
}

vector<string> FlattenFofn(const vector<string>& files)
{
    vector<string> res;

    for (const auto& file : files) {
        FlattenFofn(res, file);
    }

    return res;
}

bool ValidBaseFeatures(const PacBio::BAM::DataSet& ds)
{
    for (const auto& bam : ds.BamFiles()) {
        for (const auto& rg : bam.Header().ReadGroups()) {
            // P6-C4 and S/P1-C1/beta do not require covariates besides SNR
            if (rg.SequencingChemistry() == "P6-C4" || rg.SequencingChemistry() == "S/P1-C1/beta")
                continue;
            // everything else requires IPD and PulseWidth
            else if (!rg.HasBaseFeature(PacBio::BAM::BaseFeature::IPD) ||
                     !rg.HasBaseFeature(PacBio::BAM::BaseFeature::PULSE_WIDTH))
                return false;
        }
    }
    return true;
}

}  // namespace IO
}  // namespace PacBio
