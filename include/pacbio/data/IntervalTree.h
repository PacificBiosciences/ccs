// Author: Lance Hepler

#pragma once

#include <algorithm>
#include <set>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>

#include <pacbio/data/Interval.h>

namespace PacBio {
namespace Data {

class IntervalTree
{
private:
    struct WeakIntervalOrdering
    {
        bool operator()(const Interval& lhs, const Interval& rhs) const
        {
            return lhs.Left() < rhs.Left();
        }
    };

    typedef std::multiset<Interval, WeakIntervalOrdering> StorageType;

public:
    typedef StorageType::iterator iterator;
    typedef StorageType::const_iterator const_iterator;

    inline void Insert(const Interval& interval)
    {
        auto it = storage.insert(interval);

        // check to see if we overlap the previous element,
        // if we do, start our merge loop from there
        if (it != begin()) {
            const_iterator prev = std::prev(it);

            if (prev->Overlaps(*it)) it = prev;
        }

        while (it != end()) {
            const_iterator nx = std::next(it);

            if (nx != end() && it->Overlaps(*nx)) {
                const Interval u = it->Union(*nx);
                it = storage.erase(it);
                it = storage.erase(it);
                it = storage.insert(it, u);
            } else
                break;
        }
    }

    inline IntervalTree Gaps() const
    {
        IntervalTree gaps;

        for (auto it = begin(); it != end(); ++it) {
            const_iterator nx = std::next(it);

            if (nx == end()) break;

            gaps.Insert(Interval(it->Right(), nx->Left()));
        }

        return gaps;
    }

    inline IntervalTree Gaps(const Interval& interval) const
    {
        auto left = begin();
        auto right = end();

        // if we're empty (left == right == end()) or we don't overlap the interval
        // return just the provided interval
        if (left == right || !interval.Overlaps(Interval(left->Left(), (--right)->Right()))) {
            IntervalTree gaps;
            gaps.Insert(interval);
            return gaps;
        }

        IntervalTree gaps = Gaps();

        if (interval.Left() < left->Left()) gaps.Insert(Interval(interval.Left(), left->Left()));

        if (right->Right() < interval.Right())
            gaps.Insert(Interval(right->Right(), interval.Right()));

        return gaps;
    }

    inline bool Contains(const size_t value) const
    {
        const_iterator it = std::lower_bound(begin(), end(), Interval(value, value + 1));

        if (it != begin()) it = std::prev(it);

        for (; it != end() && it->Left() <= value; ++it) {
            if (it->Contains(value)) return true;
        }

        return false;
    }

    inline iterator begin() { return storage.begin(); }
    inline const_iterator begin() const { return storage.begin(); }
    inline iterator end() { return storage.end(); }
    inline const_iterator end() const { return storage.end(); }

    inline size_t size() const { return storage.size(); }

    static IntervalTree FromString(const std::string& str)
    {
        std::vector<std::string> components;
        boost::split(components, str, boost::is_any_of(","));
        IntervalTree tree;
        for (const auto& component : components) {
            tree.Insert(Interval::FromString(component));
        }
        return tree;
    }

private:
    StorageType storage;
};

}  // namespace Data
}  // namespace PacBio
