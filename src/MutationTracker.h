// Author: David Seifert

#pragma once

#include <algorithm>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include <boost/optional.hpp>

#include <pacbio/consensus/Mutation.h>
#include <pacbio/consensus/PolishResult.h>
#include <pacbio/data/internal/ConversionFunctions.h>

namespace PacBio {
namespace Consensus {

class MutationTracker
{
private:
    enum class TrackedMutationType : uint8_t
    {
        TEMPLATE,
        INSERTION,
        SUBSTITUTION
    };

    struct origTplInfo
    {
        int64_t origPos;
        // mutType can be either a
        //   - TEMPLATE : i.e. the original template (unchanged)
        //   - INSERTION
        //   - SUBSTITUTION : i.e. newTplBase != originalTpl[origPos]
        //
        // Notice that there is not "Deletion" type, as there is no way
        // to point from anywhere in the current Template to the original
        // one if a deletion occurred. Nonetheless, finding deletions is
        // easy, as the vector will have a discontinuity, e.g.
        //
        //            012
        //   origTpl: AAC
        //    curTpl: A-C
        //
        // where the vector will contain {{0, TEMPLATE, 'A'}, {2, TEMPLATE, 'C'}},
        // i.e., for index i and i+1, if vec[i+1].origPos - vec[i].origPos > 1
        // holds, we have lost a base.
        TrackedMutationType mutType;
        char newTplBase;
        boost::optional<double> pvalue;
    };

public:
    MutationTracker(std::string originalTpl)
        : MutationsApplied_{0}, originalTpl_{std::move(originalTpl)}
    {
        curTplToOrigTpl_.reserve(originalTpl_.size());

        const int64_t tplLength = originalTpl_.size();
        for (int64_t i = 0; i < tplLength; ++i) {
            curTplToOrigTpl_.push_back(
                {i, TrackedMutationType::TEMPLATE, originalTpl_[i], boost::none});
        }
    }

    // have to be sorted according to Mutation::SiteComparer

    // update diploid bookkeeping
    // TODO(dseifert):
    // Try and do this more implicitly with less vector rewriting and not in O(L*N).
    //
    // While in theory a linked list should be preferred for this, that is CS theory
    // and has been invalidated in practice by caches, prefetchers and branch prediction.
    // https://view.officeapps.live.com/op/view.aspx?src=http%3a%2f%2fvideo.ch9.ms%2fsessions%2fbuild%2f2014%2f2-661.pptx
    inline void AddSortedMutations(const std::vector<Mutation>& muts)
    {
        for (auto it = muts.crbegin(); it != muts.crend(); ++it) {
            // Caveat: Current diploid handling does not
            // handle Mutations having Length() > 1.
            auto insertIter = curTplToOrigTpl_.begin() + it->Start();
            int64_t tempOrigPos;
            auto tempPvalue = boost::optional<double>{};

            switch (it->Type()) {
                case MutationType::DELETION:
                    curTplToOrigTpl_.erase(curTplToOrigTpl_.begin() + it->Start(),
                                           curTplToOrigTpl_.begin() + it->Start() + it->Length());
                    break;

                case MutationType::INSERTION:
                    // save original Pos and pvalue due to
                    // inserting shifting around stuff
                    tempOrigPos = curTplToOrigTpl_[it->Start()].origPos;
                    tempPvalue = it->GetPvalue();

                    insertIter = curTplToOrigTpl_.insert(curTplToOrigTpl_.begin() + it->Start(),
                                                         it->LengthDiff(), origTplInfo{});

                    for (auto basesIter = it->Bases().cbegin(); basesIter < it->Bases().cend();
                         ++basesIter, ++insertIter) {
                        insertIter->origPos = tempOrigPos;
                        insertIter->mutType = TrackedMutationType::INSERTION;
                        insertIter->newTplBase = *basesIter;
                        insertIter->pvalue = tempPvalue;
                    }
                    break;

                case MutationType::SUBSTITUTION:
                    for (auto basesIter = it->Bases().cbegin(); basesIter < it->Bases().cend();
                         ++basesIter, ++insertIter) {
                        // only change vanilla TEMPLATE positions into SUBSTITUTIONs
                        // if a position is already an INSERTION, we need to keep
                        // tracking of the INSERTION, hence it has to stay.
                        if (insertIter->mutType == TrackedMutationType::TEMPLATE)
                            insertIter->mutType = TrackedMutationType::SUBSTITUTION;

                        insertIter->newTplBase = *basesIter;
                        insertIter->pvalue = it->GetPvalue();
                    }
                    break;
            }

            ++MutationsApplied_;
        }
    }

    // extract the meaty parts from curTplToOrigTpl_
    inline std::vector<DiploidSite> MappingToOriginalTpl() const
    {
        std::vector<DiploidSite> result;

        // leave some buffer
        result.reserve(2 * MutationsApplied_);

        // 1. find all SUBSTITUTIONs and INSERTIONs
        for (const auto& j : curTplToOrigTpl_) {
            if (j.mutType != TrackedMutationType::TEMPLATE) {
                result.emplace_back(static_cast<MutationType>(j.mutType),
                                    Data::detail::demultiplexAmbiguousBase(j.newTplBase), j.origPos,
                                    j.pvalue);
            }
        }

        // 2. find all DELETIONs
        // Recall that curTplToOrigTpl_ cannot represent deletions
        // so we find deletions by looking at all TEMPLATE positions
        // and checking whether there is a discontinuity between them
        auto lastTplPos =
            std::find_if(curTplToOrigTpl_.cbegin(), curTplToOrigTpl_.cend(), [](const auto& elem) {
                return (elem.mutType == TrackedMutationType::TEMPLATE ||
                        elem.mutType == TrackedMutationType::SUBSTITUTION);
            });
        if (lastTplPos == curTplToOrigTpl_.cend())
            throw std::runtime_error(
                "The template has been completely mutated, this should not occur!");

        // 2.1 check for deletions before the first remaining original template base
        for (int64_t delPos = 0; delPos < lastTplPos->origPos; ++delPos) {
            // deletions don't have an associated p-value
            result.emplace_back(MutationType::DELETION, std::vector<char>{}, delPos);
        }

        // 2.2 check for deletions between the first and last remaining original template bases
        for (auto it = curTplToOrigTpl_.cbegin(); it < curTplToOrigTpl_.cend(); ++it) {
            if ((it->mutType == TrackedMutationType::TEMPLATE) ||
                (it->mutType == TrackedMutationType::SUBSTITUTION)) {
                for (int64_t delPos = lastTplPos->origPos + 1; delPos < it->origPos; ++delPos) {
                    // deletions don't have an associated p-value
                    result.emplace_back(MutationType::DELETION, std::vector<char>{}, delPos);
                }
                lastTplPos = it;
            }
        }

        // 2.3 check for deletions after the last remaining original template base
        for (size_t delPos = lastTplPos->origPos + 1; delPos < originalTpl_.size(); ++delPos) {
            // deletions don't have an associated p-value
            result.emplace_back(MutationType::DELETION, std::vector<char>{},
                                static_cast<int64_t>(delPos));
        }

        // 3. finally sort everything
        std::sort(result.begin(), result.end(), [](const auto& lhs, const auto& rhs) {
            // copied from Mutation.h
            const auto l = std::make_tuple(lhs.pos, lhs.mutType != MutationType::DELETION);
            const auto r = std::make_tuple(rhs.pos, rhs.mutType != MutationType::DELETION);
            return l < r;
        });

        return result;
    }

private:
    int32_t MutationsApplied_;
    std::string originalTpl_;

    // diploid bookkeeping vector
    //
    // In order to generate the correct std::vector<Mutation>
    // for the diploid result, we need to keep track of the
    // correspondence between the current template and the
    // original one
    // This vector is a map
    //
    //   f : currentTplIdx -> [originatingIdx, +type etc]
    //
    // Thus, the length of the vector is always equal to the
    // current template.
    std::vector<origTplInfo> curTplToOrigTpl_;
};

}  // namespace Consensus
}  // namespace PacBio
