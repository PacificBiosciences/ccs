%{
/* Includes the header in the wrapper code */
#include <ConsensusCore/Sequence.hpp>
#include <ConsensusCore/Mutation.hpp>
#include <ConsensusCore/Read.hpp>
#include <ConsensusCore/Quiver/MultiReadMutationScorer.hpp>
#include <ConsensusCore/Quiver/MutationScorer.hpp>
#include <ConsensusCore/Quiver/QuiverConfig.hpp>
#include <ConsensusCore/Quiver/SimpleRecursor.hpp>
#include <ConsensusCore/Quiver/SseRecursor.hpp>
#include <ConsensusCore/Quiver/ReadScorer.hpp>
#include <ConsensusCore/Quiver/Diploid.hpp>
#include <ConsensusCore/Quiver/QuiverConsensus.hpp>

using namespace ConsensusCore;
%}

#ifdef SWIGPYTHON

%include "numpy.i"
%numpy_typemaps(float, NPY_FLOAT, int)

%apply (float* IN_ARRAY2, int DIM1, int DIM2)
       { (const float *siteScores, int dim1, int dim2) }

#endif // SWIGPYTHON

 // SWIG now seems to be incorrectly deciding that MultiReadMutationScorer
 // is an abstract class, so we have to tell it otherwise
%feature("notabstract") MultiReadMutationScorer;

#ifdef SWIGCSHARP
%csmethodmodifiers *::ToString() const "public override"
#endif // SWIGCSHARP


%include <ConsensusCore/Sequence.hpp>
%include <ConsensusCore/Mutation.hpp>
%include <ConsensusCore/Read.hpp>
%include <ConsensusCore/Quiver/detail/Combiner.hpp>
%include <ConsensusCore/Quiver/detail/RecursorBase.hpp>
%include <ConsensusCore/Quiver/MultiReadMutationScorer.hpp>
%include <ConsensusCore/Quiver/MutationScorer.hpp>
%include <ConsensusCore/Quiver/QuiverConfig.hpp>
%include <ConsensusCore/Quiver/SimpleRecursor.hpp>
%include <ConsensusCore/Quiver/SseRecursor.hpp>
%include <ConsensusCore/Quiver/ReadScorer.hpp>
%include <ConsensusCore/Quiver/Diploid.hpp>
%include <ConsensusCore/Quiver/QuiverConsensus.hpp>

 
namespace ConsensusCore {
    //
    // Dense matrix recursors and such
    //
    %template(QvRecursorBase)           detail::RecursorBase<DenseMatrix, QvEvaluator, detail::ViterbiCombiner>;
    %template(SimpleQvRecursor)         SimpleRecursor<DenseMatrix, QvEvaluator, detail::ViterbiCombiner>;
    %template(SimpleQvMutationScorer)   MutationScorer<SimpleQvRecursor>;
    %template(SseQvRecursor)            SseRecursor<DenseMatrix, QvEvaluator, detail::ViterbiCombiner>;
    %template(SseQvMutationScorer)      MutationScorer<SseQvRecursor>;

    //
    // Sparse matrix support
    //
    %template(SparseQvRecursorBase)           detail::RecursorBase<SparseMatrix, QvEvaluator, detail::ViterbiCombiner>;
    %template(SparseSimpleQvRecursor)         SimpleRecursor<SparseMatrix, QvEvaluator, detail::ViterbiCombiner>;
    %template(SparseSimpleQvMutationScorer)   MutationScorer<SparseSimpleQvRecursor>;
    %template(SparseSseQvRecursor)            SseRecursor<SparseMatrix, QvEvaluator, detail::ViterbiCombiner>;
    %template(SparseSseQvMutationScorer)      MutationScorer<SparseSseQvRecursor>;

    %template(SparseSseQvMultiReadMutationScorer) MultiReadMutationScorer<SparseSseQvRecursor>;

    //
    // Sparse matrix sum-product support
    //
    %template(SparseQvSumProductRecursorBase)           detail::RecursorBase<SparseMatrix, QvEvaluator, detail::SumProductCombiner>;
    %template(SparseSimpleQvSumProductRecursor)         SimpleRecursor<SparseMatrix, QvEvaluator, detail::SumProductCombiner>;
    %template(SparseSimpleQvSumProductMutationScorer)   MutationScorer<SparseSimpleQvSumProductRecursor>;
    %template(SparseSseQvSumProductRecursor)            SseRecursor<SparseMatrix, QvEvaluator, detail::SumProductCombiner>;
    %template(SparseSseQvSumProductMutationScorer)      MutationScorer<SparseSseQvSumProductRecursor>;

    %template(SparseSseQvSumProductMultiReadMutationScorer) MultiReadMutationScorer<SparseSseQvSumProductRecursor>;

    //
    // Edna evaluator support
    //
    %template(SparseEdnaRecursorBase)           detail::RecursorBase<SparseMatrix, EdnaEvaluator, detail::SumProductCombiner>;
    %template(SparseSseEdnaRecursor)            SseRecursor<SparseMatrix, EdnaEvaluator, detail::SumProductCombiner>;
    %template(SparseSseEdnaMutationScorer)      MutationScorer<SparseSseEdnaRecursor>;
}
