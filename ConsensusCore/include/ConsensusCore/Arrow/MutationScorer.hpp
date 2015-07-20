
#pragma once

#include <boost/noncopyable.hpp>
#include <string>

// TODO(dalexander): how can we remove this include??
//  We should move all template instantiations out to another
//  header, I presume.
#include <ConsensusCore/Mutation.hpp>
#include <ConsensusCore/Types.hpp>
#include <ConsensusCore/Arrow/SimpleRecursor.hpp>
#include <ConsensusCore/Arrow/ContextParameters.hpp>

namespace ConsensusCore {
namespace Arrow {
   
    template<typename R>
    class MutationScorer
    {
        public:
            typedef typename R::MatrixType    MatrixType;
            typedef R                         RecursorType;

        public:
            MutationScorer(const R& recursor)
                throw(AlphaBetaMismatchException);

            MutationScorer(const MutationScorer& other);
            virtual ~MutationScorer();

        public:
            WrappedTemplateParameterPair Template() const;
            void Template(WrappedTemplateParameterPair tpl)
                throw(AlphaBetaMismatchException);

            double Score() const;
            // Must be called from a multiread mutation context.
            double ScoreMutation(const Mutation& m) const;

        private:
            void DumpMatrix(const MatrixType& mat, const std::string& fname) const;

        public:
            void DumpBetaMatrix() const;
            void DumpAlphaMatrix() const;

        public:
            // Accessors that are handy for debugging.
            const MatrixType* Alpha() const;
            const MatrixType* Beta() const;
            const int NumFlipFlops() const { return numFlipFlops_; }

        private:
             R* recursor_;
            /**
             The Alpha matrix, holding the forward variables
             */
            MatrixType* alpha_;
            MatrixType* beta_;
            /**
             The extension matrix, a clean buffer to work with data where we need it.
             */
            MatrixType* extendBuffer_;
            int numFlipFlops_;
        
        
            //friend class MultiReadMutationScorer<R>;
    };
    
    typedef MutationScorer<ArrowRecursor> ArrowMutationScorer;
}    
}
