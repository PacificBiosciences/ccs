// Author: Lance Hepler

#include <pacbio/consensus/IntervalMask.h>

namespace PacBio {
namespace Consensus {
namespace {

size_t SafeAdd(size_t a, int b)
{
    if (b < 0 && static_cast<size_t>(-b) >= a) return static_cast<size_t>(0);
    return a + b;
}
}

bool IntervalMask::Contains(const Mutation& mut)
{
    if (mut.Type() == MutationType::INSERTION)
        return IntervalTree::Contains(mut.End()) &&
               (mut.End() == 0 || IntervalTree::Contains(mut.End() - 1));
    return IntervalTree::Contains(mut.Start());
}

void IntervalMask::Mutate(const std::vector<Mutation>& muts)
{
    if (muts.empty()) return;
    // iterate through each interval in the mask, and:
    //   1) deplete mutations (if any) to the left of the interval, updating offL
    //   2) deplete mutations (if any) within the interval, updating offR
    //   3) add the interval to the new mask, upating (L, R) w/ offL and offR respectively
    IntervalMask newMask;
    auto m = muts.begin();
    int offL = 0;
    for (auto ab : *this) {
        //   if the mutation's right is before the interval's left, update offL
        for (; m != muts.end() && m->End() <= ab.Left(); ++m)
            offL += m->LengthDiff();
        int offR = offL;
        // if the mutation's left is within the interval, then update offR
        for (; m != muts.end() && ab.Contains(m->Start()); ++m)
            offR += m->LengthDiff();
        size_t l = SafeAdd(ab.Left(), offL);
        size_t r = SafeAdd(ab.Right(), offR);
        // if the interval has a span, add it
        if (l < r) newMask.Insert({l, r});
        // offR is the new offL
        offL = offR;
    }
    (*this) = std::move(newMask);
}
}
}
