
#include <algorithm>
#include <numeric>
#include <cmath>

#include <boost/tuple/tuple.hpp>

#include <ConsensusCore/Interval.hpp>
#include <ConsensusCore/LValue.hpp>


namespace ConsensusCore {

template<typename M>
ScaledMatrix<M>::ScaledMatrix(int rows, int cols)
    : M(rows, cols)
    , logScalars_(cols, static_cast<F>(0))
{ }

template<typename M>
ScaledMatrix<M>::ScaledMatrix(const ScaledMatrix<M>& other)
    : M(other)
    , logScalars_(other.logScalars_)
{ }

template<typename M>
inline const ScaledMatrix<M>&
ScaledMatrix<M>::Null()
{
    static ScaledMatrix<M>* nullObj = new ScaledMatrix<M>(0, 0);
    return *nullObj;
}

template<typename M>
inline void
ScaledMatrix<M>::FinishEditingColumn(int j, int usedBegin, int usedEnd)
{
    // get the constant to scale by
    F c = static_cast<F>(0);
    for (int i = usedBegin; i < usedEnd; ++i)
    {
        c = std::max(c, M::Get(i, j));
    }

    // set it
    if (c != static_cast<F>(0) && c != static_cast<F>(1))
    {
        for (int i = usedBegin; i < usedEnd; ++i)
        {
            M::Set(i, j, M::Get(i, j) / c);
        }
        logScalars_[j] = std::log(c);
    }
    else
    {
        logScalars_[j] = static_cast<F>(0);
    }

    M::FinishEditingColumn(j, usedBegin, usedEnd);
}

template<typename M>
inline typename M::FloatType
ScaledMatrix<M>::GetLogScale(int j) const
{
    return logScalars_[j];
}

template<typename M>
inline typename M::FloatType
ScaledMatrix<M>::GetLogProdScales(int beginColumn, int endColumn) const
{
    return std::accumulate(logScalars_.begin() + beginColumn,
                           logScalars_.begin() + endColumn,
                           static_cast<F>(0));
}

template<typename M>
inline typename M::FloatType
ScaledMatrix<M>::GetLogProdScales() const
{
    return std::accumulate(logScalars_.begin(),
                           logScalars_.end(),
                           static_cast<F>(0));
}

}
