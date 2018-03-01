// Author: Brett Bowman

#pragma once

#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <memory>
#include <numeric>

// Initialize data structures, do NOT remove
#include <pacbio/consensus/internal/ModelInternalInitializer.h>

#include <pacbio/consensus/Evaluator.h>
#include <pacbio/consensus/Mutation.h>
#include <pacbio/data/Read.h>
#include <pacbio/data/State.h>
#include <pacbio/exception/StateError.h>

namespace PacBio {
namespace Consensus {

// Forward decl
class AbstractMatrix;

/// Contains user-provided filtering information for the Evaluators.
struct IntegratorConfig
{
    double MinZScore;
    double ScoreDiff;

    IntegratorConfig(double minZScore = -3.4, double scoreDiff = 25.0);
};

/// The Integrator holds a collection of Evaluators whose MappedReads belonging
/// to the same genomic region or amplicon
class Integrator
{
public:
    /// \brief Initialize the Integrator.
    ///
    /// \param tpl    The draft template as a string
    /// \param cfg    The configuration used to initialize the Integrator.
    Integrator(const std::string& tpl, const IntegratorConfig& cfg);

    Integrator(Integrator&&) = default;

public:
    virtual ~Integrator() = default;

    virtual size_t TemplateLength() const;

    /// Returns base i of the template
    virtual char operator[](size_t i) const;
    virtual operator std::string() const;

    /// This method throws InvalidEvaluatorException, every time the likelihood
    /// can't be computed for one Evaluator; the respective Evaluator will be invalidated.
    /// You MUST recompute the LLs for all your mutations of interest, as the
    /// number of active Evaluators changed.
    virtual double LL(const Mutation& mut);
    virtual double LL() const;

    /// Masks intervals of the template for each read where the observed error rate is
    /// greater than maxErrRate in 1+2*radius template bases
    void MaskIntervals(size_t radius, double maxErrRate);

    /// Applies a mutation to the template of each Evaluator.
    virtual void ApplyMutation(const Mutation& mut);
    /// Applies a vector of murations to the template of each Evaluator.
    virtual void ApplyMutations(std::vector<Mutation>* muts);

    /// Encapsulate the read in an Evaluator and stores it.
    virtual PacBio::Data::State AddRead(const PacBio::Data::MappedRead& read);

public:
    double AvgZScore() const;
    std::vector<double> ZScores() const;
    std::vector<std::pair<double, double>> NormalParameters() const;

    /// Given a Mutation of interest, returns a vector of LLs,
    /// one LL per active Evaluator; invalid Evaluators are omitted.
    ///
    /// This method throws InvalidEvaluatorException, every time the likelihood
    /// can't be computed for one Evaluator; the respective Evaluator will be invalidated.
    /// You MUST recompute the LLs for all your mutations of interest, as the
    /// number of active Evaluators changed.
    std::vector<double> LLs(const Mutation& mut);
    /// Return the LL for each Evaluator, even invalid ones.
    /// DO NOT use this in production code, only for debugging purposes.
    std::vector<double> LLs() const;
    /// Return the best-mutation improvement histogram for a locus
    /// and given MutationType.
    ///
    /// Say we have 10 Evaluators, and provide some site and a MutationType::INSERTION
    ///   - 3 of them, 'A' provides the best LL improvement,
    ///   - 0 of them, 'C' provides the best LL improvement,
    ///   - 1 of them, 'G' provides the best LL improvement,
    ///   - 5 of them, 'T' provides the best LL improvement.
    /// Notice that the sum of these is 9, that is, 1 Evaluator is either invalid
    /// or its LL decreases for every base.
    /// The return value is a reverse-sorted array of the base and the number of Evaluators,
    /// sorted on the second field, e.g., for the aforementioned example we will have
    ///
    ///   {{'T', 5}, {'A', 3}, {'G', 1}, {'C', 0}}
    ///
    /// as return value.
    std::array<std::pair<char, int>, 4> BestMutationHistogram(size_t start, MutationType mutType);
    /// For each Evaluator, returns the read name.
    std::vector<std::string> ReadNames() const;
    /// Returns the number of flip flop events for each Evaluator.
    std::vector<int> NumFlipFlops() const;
    /// Returns the maximal number of flip flop events of all Evaluators.
    int MaxNumFlipFlops() const;
    /// Computes the ratio of populated cells in the alpha matrix for each
    /// Evaluator and returns the maximal ratio.
    float MaxAlphaPopulated() const;
    /// Computes the ratio of populated cells in the beta matrix for each
    /// Evaluator and returns the maximal ratio.
    float MaxBetaPopulated() const;
    /// Returns the state of each Evaluator.
    std::vector<PacBio::Data::State> States() const;
    /// Returns the strand of each Evaluator.
    std::vector<PacBio::Data::StrandType> StrandTypes() const;

    /// Returns read-only access to Evaluator idx.
    const Evaluator& GetEvaluator(size_t idx) const;

public:
    // Abstract matrix access for SWIG and diagnostics
    const AbstractMatrix& Alpha(size_t idx) const;
    const AbstractMatrix& Beta(size_t idx) const;

protected:
    Mutation ReverseComplement(const Mutation& mut) const;

    PacBio::Data::State AddRead(std::unique_ptr<AbstractTemplate>&& tpl,
                                const PacBio::Data::MappedRead& read);

    std::unique_ptr<AbstractTemplate> GetTemplate(const PacBio::Data::MappedRead& read);

protected:
    IntegratorConfig cfg_;
    std::vector<Evaluator> evals_;
    std::string fwdTpl_;
    std::string revTpl_;

private:
    /// Return LL for a single Evaluator
    template <bool AllowInvalidEvaluators>
    inline double SingleEvaluatorLL(Evaluator* const eval, const Mutation& fwdMut) const;

    /// Extract a feature vector from a vector of Evaluators for non-const functions.
    template <typename T>
    inline std::vector<T> TransformEvaluators(std::function<T(Evaluator&)> functor)
    {
        // TODO(atoepfer) How can we use const_cast to convert this?
        std::vector<T> vec;
        vec.reserve(evals_.size());
        std::transform(evals_.begin(), evals_.end(), std::back_inserter(vec), functor);
        return vec;
    }

    /// Extract a feature vector from a vector of Evaluators for const functions.
    template <typename T>
    inline std::vector<T> TransformEvaluators(std::function<T(const Evaluator&)> functor) const
    {
        std::vector<T> vec;
        vec.reserve(evals_.size());
        std::transform(evals_.begin(), evals_.end(), std::back_inserter(vec), functor);
        return vec;
    }

private:
    friend struct std::hash<Integrator>;
};

/// Helper function to get maximal number from a vector.
template <typename T>
inline T MaxElement(const std::vector<T>& in)
{
    return *std::max_element(in.cbegin(), in.end());
}

}  // namespace Consensus
}  // namespace PacBio
