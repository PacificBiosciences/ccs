// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

#include <pacbio/ccs/SubreadResultCounter.hpp>
#include <stdexcept>

using namespace PacBio::Consensus;
using namespace PacBio::CCS;

std::vector<int32_t> SubreadResultCounter::ReturnCountsAsArray() const {
        std::vector<int32_t> results {Success, AlphaBetaMismatch, BelowMinQual, PoorZScore,
            FilteredBySize, Other};
    return results;
}

void SubreadResultCounter::AddResult(AddReadResult result) {
    switch (result) {
        case AddReadResult::ALPHA_BETA_MISMATCH :
            AlphaBetaMismatch++;
            break;
        case AddReadResult::OTHER:
            Other++;
            break;
        case AddReadResult::POOR_ZSCORE:
            PoorZScore++;
            break;
        case AddReadResult::SIZE:
            FilteredBySize++;
            break;
        case AddReadResult::SUCCESS:
            Success++;
            break;
        default:
            throw std::runtime_error("Unexpected AddReadResult");
            break;
    }
}

int32_t SubreadResultCounter::Total() const {
    return (AlphaBetaMismatch + Success + BelowMinQual + FilteredBySize + Other + PoorZScore + ZMWBelowMinSNR + ZMWNotEnoughSubReads);
}

SubreadResultCounter::SubreadResultCounter() :
        Success{0},
        AlphaBetaMismatch{0},
        BelowMinQual{0},
        FilteredBySize{0},
        ZMWBelowMinSNR{0},
        ZMWNotEnoughSubReads{0},
        PoorZScore{0},
        Other{0}

{}


void SubreadResultCounter::WriteResultsReport(std::ostream& report) const
{
    using namespace std;
    double total = static_cast<float>(Total());    

    report << "Subread Yield" << endl;
    
    report << "Success - Used for CCS," <<  Success << "," << 100.0 * Success / total
    << '%' << endl;
    
    report << "Failed -- Below SNR threshold," << ZMWBelowMinSNR << ","
    << 100.0 * ZMWBelowMinSNR / total << '%' << endl;
    
    report << "Failed -- Alpha/Beta mismatch," << AlphaBetaMismatch << ","
    << 100.0 * AlphaBetaMismatch / total << '%' << endl;
    
    report << "Failed -- Below minimum quality," << BelowMinQual << ","
    << 100.0 * BelowMinQual / total << '%' << endl;
    
    report << "Failed -- Filtered by size," << FilteredBySize << ","
    << 100.0 * FilteredBySize / total << '%' << endl;
    
    report << "Failed -- Z-Score too low," << PoorZScore << ","
    << 100.0 * PoorZScore / total << '%' << endl;
    
    report << "Failed -- From ZMW with too few passes," << ZMWNotEnoughSubReads << ","
    << 100.0 * ZMWNotEnoughSubReads / total << '%' << endl;
    
    report << "Failed -- Other," << Other << ","
    << 100.0 * Other / total << '%' << endl;
}


void SubreadResultCounter::AssignSuccessToOther() {
    Other += Success;
    Success = 0;
}

void SubreadResultCounter::CombineWithOtherResult(const SubreadResultCounter& other) {
    Success += other.Success;
    AlphaBetaMismatch += other.AlphaBetaMismatch;
    BelowMinQual += other.BelowMinQual;
    FilteredBySize += other.FilteredBySize;
    Other += other.Other;
    PoorZScore += other.PoorZScore;
    ZMWBelowMinSNR += other.ZMWBelowMinSNR;
    ZMWNotEnoughSubReads += other.ZMWNotEnoughSubReads;
}

SubreadResultCounter& SubreadResultCounter::operator+=(const SubreadResultCounter& other)
{
    CombineWithOtherResult(other);
    return *this;
}

