// Author: Nigel Delaney

#include <pacbio/data/SubreadResultCounter.h>
#include <stdexcept>

using std::string;
using std::vector;

using std::endl;

using namespace PacBio::Data;

std::vector<int32_t> SubreadResultCounter::ReturnCountsAsArray() const
{
    std::vector<int32_t> results{Success,    AlphaBetaMismatch, BelowMinQual,
                                 PoorZScore, FilteredBySize,    Other};
    return results;
}

void SubreadResultCounter::AddResult(State result)
{
    switch (result) {
        case State::ALPHA_BETA_MISMATCH:
            AlphaBetaMismatch++;
            break;
        case State::POOR_ZSCORE:
            PoorZScore++;
            break;
        case State::VALID:
            Success++;
            break;
        default:
            Other++;
            break;
    }
}

int32_t SubreadResultCounter::Total() const
{
    return (AlphaBetaMismatch + Success + BelowMinQual + FilteredBySize + Other + PoorZScore +
            ZMWBelowMinSNR + ZMWNotEnoughSubReads);
}

SubreadResultCounter::SubreadResultCounter()
    : Success{0}
    , AlphaBetaMismatch{0}
    , BelowMinQual{0}
    , FilteredBySize{0}
    , ZMWBelowMinSNR{0}
    , ZMWNotEnoughSubReads{0}
    , PoorIdentity{0}
    , PoorZScore{0}
    , Other{0}

{
}

void SubreadResultCounter::WriteResultsReport(std::ostream& report) const
{
    double total = static_cast<float>(Total());

    report << "Subread Yield" << endl;

    report << "Success - Used for CCS," << Success << "," << 100.0 * Success / total << '%' << endl;

    report << "Failed -- Below SNR threshold," << ZMWBelowMinSNR << ","
           << 100.0 * ZMWBelowMinSNR / total << '%' << endl;

    report << "Failed -- Alpha/Beta mismatch," << AlphaBetaMismatch << ","
           << 100.0 * AlphaBetaMismatch / total << '%' << endl;

    report << "Failed -- Below minimum quality," << BelowMinQual << ","
           << 100.0 * BelowMinQual / total << '%' << endl;

    report << "Failed -- Filtered by size," << FilteredBySize << ","
           << 100.0 * FilteredBySize / total << '%' << endl;

    report << "Failed -- Identity too low," << PoorIdentity << "," << 100.0 * PoorIdentity / total
           << '%' << endl;

    report << "Failed -- Z-Score too low," << PoorZScore << "," << 100.0 * PoorZScore / total << '%'
           << endl;

    report << "Failed -- From ZMW with too few passes," << ZMWNotEnoughSubReads << ","
           << 100.0 * ZMWNotEnoughSubReads / total << '%' << endl;

    report << "Failed -- Other," << Other << "," << 100.0 * Other / total << '%' << endl;
}

void SubreadResultCounter::AssignSuccessToOther()
{
    Other += Success;
    Success = 0;
}

void SubreadResultCounter::CombineWithOtherResult(const SubreadResultCounter& other)
{
    Success += other.Success;
    AlphaBetaMismatch += other.AlphaBetaMismatch;
    BelowMinQual += other.BelowMinQual;
    FilteredBySize += other.FilteredBySize;
    Other += other.Other;
    PoorIdentity += other.PoorIdentity;
    PoorZScore += other.PoorZScore;
    ZMWBelowMinSNR += other.ZMWBelowMinSNR;
    ZMWNotEnoughSubReads += other.ZMWNotEnoughSubReads;
}

SubreadResultCounter& SubreadResultCounter::operator+=(const SubreadResultCounter& other)
{
    CombineWithOtherResult(other);
    return *this;
}
