// Author: Lance Hepler

#pragma once

#include <algorithm>
#include <cassert>
#include <vector>

namespace PacBio {
namespace Poa {
namespace detail {

// This is necessary because of dumb C++ linker rules.
// see: http://isocpp.org/wiki/faq/templates#template-friends
// and: http://yosefk.com/c++fqa/templates.html#fqa-35.16  (for hilarious
// rejoinder)
template <typename T>
class VectorL;  // fwd
template <typename T>
T Max(const VectorL<T>& v);
template <typename T>
size_t ArgMax(const VectorL<T>& v);

//
// Vector class that stores only a subsequence of the rows
// (SparseVector has a lot of weird functionality built-in, so we can't use it
// without
//  cleaning it up/refactoring it quite a bit)
//
template <typename T>
class VectorL
{
private:
    std::vector<T> storage_;
    size_t beginRow_;
    size_t endRow_;

public:
    VectorL(int beginRow, int endRow, T defaultVal = T())
        : storage_(endRow - beginRow, defaultVal), beginRow_(beginRow), endRow_(endRow)
    {
    }

    T& operator[](size_t pos)
    {
        assert(beginRow_ <= pos && pos < endRow_);
        return storage_[pos - beginRow_];
    }

    const T& operator[](size_t pos) const
    {
        assert(beginRow_ <= pos && pos < endRow_);
        return storage_[pos - beginRow_];
    }

    size_t BeginRow() const { return beginRow_; }
    size_t EndRow() const { return endRow_; }
    friend T Max<>(const VectorL<T>& v);
    friend size_t ArgMax<>(const VectorL<T>& v);
};

using std::max_element;
using std::distance;

template <typename T>
T Max(const VectorL<T>& v)
{
    return *max_element(v.storage_.begin(), v.storage_.end());
}

template <typename T>
size_t ArgMax(const VectorL<T>& v)
{
    return v.beginRow_ +
           distance(v.storage_.begin(), max_element(v.storage_.begin(), v.storage_.end()));
}

}  // namespace detail
}  // namespace Poa
}  // PacBio
