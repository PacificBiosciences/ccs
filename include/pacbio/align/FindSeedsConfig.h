// Author: Lance Hepler

#pragma once

#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/sequence.h>

namespace PacBio {
namespace Align {

/// A simple struct for represting a complete specialization for a SeqAn index
/// such that template-ifying code around different indices can be done around
/// a single, simple parameter.
///
/// TODO (lhepler) : investigate default values other than 10
///
/// \param  size_t  The length of the K-mer or Q-Gram used by the index.
///                 By default 10.
/// \param  TShape  The SeqAn shape specialization.  By default Ungapped.
/// \param  TIndex  The SeqAn index specialization.  By default Q-Gram.
template <size_t TSize = 10, typename TShape = seqan::UngappedShape<TSize>,
          typename TIndex = seqan::IndexQGram<TShape>>
struct FindSeedsConfig
{
    typedef TIndex IndexType;
    typedef TShape ShapeType;
    static const size_t Size = TSize;
};

}  // Align
}  // PacBio
