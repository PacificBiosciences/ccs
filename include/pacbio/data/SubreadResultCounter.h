// Author: Armin Töpfer

#pragma once

#include <stdlib.h>
#include <iostream>
#include <vector>

#include <pacbio/data/State.h>

namespace PacBio {
namespace Data {

// A class to store and report on the results of what happens to subreads.
class SubreadResultCounter
{
public:
    int32_t Success;
    int32_t AlphaBetaMismatch;
    int32_t BelowMinQual;
    int32_t FilteredBySize;
    int32_t ZMWBelowMinSNR;
    int32_t ZMWNotEnoughSubReads;
    int32_t PoorIdentity;
    int32_t PoorZScore;
    int32_t Other;

    SubreadResultCounter();
    std::vector<int32_t> ReturnCountsAsArray() const;
    void AddResult(PacBio::Data::State);
    /* Certain conditions may make reads that were on their
       way to success go to the garbage bin, in this case we reassign
       all the success reads to the other category */
    void AssignSuccessToOther();
    void CombineWithOtherResult(const SubreadResultCounter& other);
    SubreadResultCounter& operator+=(const SubreadResultCounter& other);
    void WriteResultsReport(std::ostream& report) const;
    int32_t Total() const;
};

}  // namespace Data
}  // namespace PacBio
