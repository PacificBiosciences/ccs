// Copyright (c) 2011-2016, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

#pragma once

#include <memory>
#include <utility>
#include <vector>

#include <pacbio/consensus/MatrixViewConvention.h>
#include <pacbio/consensus/Template.h>
#include <pacbio/data/Read.h>
#include <pacbio/data/State.h>

namespace PacBio {
namespace Consensus {

// forward declaration
class EvaluatorImpl;
class AbstractMatrix;

/// Each Evaluator holds one reference to a MappedRead and its Template.
/// An Evaluator can compute the LL that its MappedRead stems from the Template.
/// Core functionality: compute the LL given a temporary mutated Template
/// and apply mutations to the Template.
///
/// A PIMPL wrapper around the implementation of the Evaluator allows to
/// deactivate the instance, either implicitly by an error or
/// expliciltly by releasing the implementation pointer.
///
/// If a function gets called on an deactivated Evaluator, it returns -INF!
class Evaluator
{
public:
    Evaluator() = delete;

    /// Initializes an empty instance as a placeholder.
    explicit Evaluator(PacBio::Data::State);

    /// Default constructor.
    ///
    /// \param tpl        The respective template.
    /// \param mr         The MappedRead
    /// \param minZScore  The minimum z-score
    /// \param scoreDiff  The score difference
    Evaluator(std::unique_ptr<AbstractTemplate>&& tpl, const PacBio::Data::MappedRead& mr,
              double minZScore, double scoreDiff);

    /// Copying is verboten
    Evaluator(const Evaluator&) = delete;
    Evaluator& operator=(const Evaluator&) = delete;

    /// Move constructor
    Evaluator(Evaluator&&);
    /// Move assign operator
    Evaluator& operator=(Evaluator&&);

    /// Destructor
    ~Evaluator();

    size_t Length() const;  // TODO: is this used anywhere?  If not, delete it.

    /// Returns the strand if Evaluator.
    /// Return StrandType::UNMAPPED if deactivated.
    PacBio::Data::StrandType Strand() const;

    /// Returns if the Evaluator is valid.
    operator bool() const { return IsValid(); }

    /// Returns if the Evaluator is still valid and active.
    bool IsValid() const { return curState_ == PacBio::Data::State::VALID; }

    /// TODO: Attention, not implemented, thus can be removed.
    operator std::string() const;

    /// Returns the read name.
    /// Returns *Inactive evaluator* if deactivated.
    std::string ReadName() const;

    /// Returns the LL of the Read, given the mutated template.
    /// Returns -INF if deactivated.
    double LL(const Mutation& mut);

    /// Returns the LL of the Read, given the current template.
    /// Returns -INF if deactivated.
    double LL() const;

    /// Returns the mean and variance over all site-wise normal parameters.
    /// Returns {-INF, -INF} if deactivated.
    std::pair<double, double> NormalParameters() const;

    /// Returns the ZScore of this Evaluator's LL, given all Evaluators of
    /// the template.
    /// Returns -INF if deactivated.
    double ZScore() const;

    /// Applies a single mutation to the template.
    /// Returns if mutation has been applied
    /// and deactivates the Evaluator if not.
    bool ApplyMutation(const Mutation& mut);

    /// Applies a vector of mutations to the template.
    /// Returns if all mutations have been applied
    /// and deactivates the Evaluator if not.
    bool ApplyMutations(std::vector<Mutation>* muts);

    /// Returns the current state of the Evaluator.
    PacBio::Data::State Status() const { return curState_; }

    /// Returns number of flip flop events from the initial alpha/beta fill.
    /// Returns -INF if deactivated.
    int NumFlipFlops() const;

    /// Invalides this Evaluator and releases its implementation.
    void Release();

public:
    const AbstractMatrix& Alpha() const;
    const AbstractMatrix& Beta() const;

    const AbstractMatrix* AlphaView(MatrixViewConvention c) const;
    const AbstractMatrix* BetaView(MatrixViewConvention c) const;

private:
    /// Checks the z-score and disables the Evaluator if does not pass
    /// the threshold.
    /// This filter noops for Sequel models.
    void CheckZScore(const double minZScore, const std::string& model);

    /// Sets the state of the Evaluator.
    /// Allows transition from VALID to anything and from anything to DISABLED.
    /// Disables the Evaluator if not VALID.
    void Status(PacBio::Data::State nextState);

private:
    std::unique_ptr<EvaluatorImpl> impl_;
    PacBio::Data::State curState_;
};

}  // namespace Consensus
}  // namespace PacBio
