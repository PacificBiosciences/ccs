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

#include <cassert>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include <pacbio/consensus/Template.h>
#include <pacbio/exception/StateError.h>

namespace PacBio {
namespace Consensus {

using TemplateTooSmall = PacBio::Exception::TemplateTooSmall;

//
// AbstractTemplate Function Definitions
//
AbstractTemplate::AbstractTemplate(const size_t start, const size_t end, const bool pinStart,
                                   const bool pinEnd)
    : start_{start}, end_{end}, pinStart_{pinStart}, pinEnd_{pinEnd}
{
    assert(start_ <= end_);

    // Templates that are below two bases are not allowed.
    // This exception will get propagated up the chain up to the Integrator
    // constructor or Integrator::AddRead.
    if (end_ - start_ < 2) throw TemplateTooSmall();
}

AbstractTemplate::~AbstractTemplate() {}
AbstractTemplate::operator std::string() const
{
    std::stringstream ss;
    for (size_t i = 0; i < Length(); ++i)
        ss << (*this)[i].Base;
    return ss.str();
}

bool AbstractTemplate::ApplyMutation(const Mutation& mut)
{
    const bool mutApplied = InRange(mut.Start(), mut.End());

    // update the end_ mapping if...
    //   we're pinned at the end (and we're not trying to delete nonexistent bases), or
    //   we're before the end of our mapping, or
    //   TODO(lhepler) I HAVE NO IDEA WHY I PUT THIS HERE ARGH
    if ((pinEnd_ && (end_ > 0 || mut.LengthDiff() > 0)) || mut.Start() < end_ ||
        mut.End() <= start_)
        end_ += mut.LengthDiff();

    // update the start_ mapping if...
    //   we're NOT pinned at the start, or
    //   the mutation comes before the start mapping.
    if (!pinStart_ && mut.End() <= start_) start_ += mut.LengthDiff();

    assert(start_ <= end_);

    return mutApplied;
}

bool AbstractTemplate::ApplyMutations(std::vector<Mutation>* const muts)
{
    bool mutsApplied = false;

    // make sure the mutations are sorted by site: End() then Start()
    std::sort(muts->begin(), muts->end(), Mutation::SiteComparer);

    for (auto it = muts->crbegin(); it != muts->crend(); ++it)
        mutsApplied |= ApplyMutation(*it);

    return mutsApplied;
}

std::pair<double, double> AbstractTemplate::NormalParameters() const
{
    double mean = 0.0, var = 0.0;
    /* Brute forcing this for clarity of code as I don't think this is limiting.
     Should profiling show this is a problem, we can use the fact that there are only
     8 possible dinucleotide contexts and so we can multiply instead of performing repeated
     addition.
     e.g. A + B + A + A + B = 3 * A + 2 * B, so we could avoid the nested add step */

    // Now sum up the mean and variance
    for (size_t i = 0; (i + 1) < Length(); ++i) {
        double m, v;
        std::tie(m, v) = SiteNormalParameters(i);
        mean += m;
        /* Add the variance, note that the sites are independent so
           Var (A + B) = Var(A) + Var(B) as Cov(A,B) = 0 */
        var += v;  // Add the variance (note sites are indpe
    }
    return std::make_pair(mean, var);
}

/* See PBEP #4 for a write-up of this code and an explanation of the algorithm.

     The R script below can be used to validate the moments are calculated correctly by
     comparing to a brute-force simulation

     # Sample Parameters
     p_m  = 0.95583140484751283
     p_d  = 0.00097238955012494488
     p_b  = 0.029256323818866534
     p_s  = 0.013939881783495679
     eps  = 0.00505052456472967

     # Expected Results
     mean  = -0.27568172991312162
     var = 1.019204780302317

     pmE = p_m / (p_m + p_d)
     exitLL <- function() {
     if (runif(1) < pmE) {
     if (runif(1) < eps) {
     return( log(p_m) + log(eps) + log(1/3))
     } else {
     return(log(p_m) + log(1-eps))
     }
     } else {
     return(log(p_d))
     }
     }

     insertLL <- function() {
     LL <- 0
     pbI = p_b / (p_b + p_s)
     while (runif(1) < (p_b + p_s)) {
     if (runif(1) < pbI) {
     LL <- LL + log(p_b)
     } else {
     LL <- LL + log(p_s) + log(1/3)
     }
     }
     return(LL)
     }

     getSamp <- function() {
     return(insertLL() + exitLL())
     }
     res = replicate(5000000, getSamp())
     mean(res)
     var(res)
*/
std::pair<double, double> AbstractTemplate::SiteNormalParameters(const size_t i) const
{
    const auto params = (*this)[i];
    // TODO: This is a bit unsafe, we should understand context and have this conversion in one
    // place (or really just move directly to using .Idx without bit shifting
    const uint8_t prev = (i == 0) ? 0 : (*this)[i - 1].Idx;  // default base : A
    const uint8_t curr = params.Idx;

    const double p_m = params.Match, l_m = std::log(p_m), l2_m = l_m * l_m;
    const double p_d = params.Deletion, l_d = std::log(p_d), l2_d = l_d * l_d;
    const double p_b = params.Branch, l_b = std::log(p_b), l2_b = l_b * l_b;
    const double p_s = params.Stick, l_s = std::log(p_s), l2_s = l_s * l_s;
    const double p_n = p_m + p_d;  // next
    const double p_e = p_b + p_s;  // extra

    // first moment expectations (zero terms used for clarity)
    const double E_M = ExpectedLLForEmission(MoveType::MATCH, prev, curr, MomentType::FIRST);
    const double E_D = 0.0;
    const double E_B = ExpectedLLForEmission(MoveType::BRANCH, prev, curr, MomentType::FIRST);
    const double E_S = ExpectedLLForEmission(MoveType::STICK, prev, curr, MomentType::FIRST);
    const double E_N = (l_m + E_M) * p_m / p_n + (l_d + E_D) * p_d / p_n;
    const double E_E = (l_b + E_B) * p_b / p_e + (l_s + E_S) * p_s / p_e;

    // calculate first moment
    const double mean = E_N + p_e * E_E / p_n;

    // second moment expectations
    const double E2_M = ExpectedLLForEmission(MoveType::MATCH, prev, curr, MomentType::SECOND);
    const double E2_D = 0.0;
    const double E2_S = ExpectedLLForEmission(MoveType::STICK, prev, curr, MomentType::SECOND);
    const double E2_B = ExpectedLLForEmission(MoveType::BRANCH, prev, curr, MomentType::SECOND);
    const double E2_N = (l2_m + 2 * l_m * E_M + E2_M) * p_m / p_n + l2_d * p_d / p_n;
    const double E2_E =
        (l2_b + 2 * E_B * l_b + E2_B) * p_b / p_e + (l2_s + 2 * E_S * l_s + E2_S) * p_s / p_e;

    // calculate second moment
    const double E2_LL = E2_N + 2 * p_e * E_N * E_E / p_n + p_e * (1 + p_e) * E2_E / (p_n * p_n);
    const double var = E2_LL - (mean * mean);

    return std::make_pair(mean, var);
}

bool AbstractTemplate::InRange(const size_t start, const size_t end) const
{
    if ((pinStart_ || start_ < end) && (pinEnd_ || start < end_)) return true;
    return false;
}

std::ostream& operator<<(std::ostream& os, const AbstractTemplate& tpl)
{
    os << std::string(tpl);
    return os;
}

//
// (Concrete) Template Function Definitions
//
Template::Template(const std::string& tpl, std::unique_ptr<ModelConfig>&& cfg)
    : Template(tpl, std::forward<std::unique_ptr<ModelConfig>>(cfg), 0, tpl.length(), true, true)
{
}

Template::Template(const std::string& tpl, std::unique_ptr<ModelConfig>&& cfg, const size_t start,
                   const size_t end, const bool pinStart, const bool pinEnd)
    : AbstractTemplate(start, end, pinStart, pinEnd)
    , cfg_(std::forward<std::unique_ptr<ModelConfig>>(cfg))
    , tpl_{cfg_->Populate(tpl)}
{
    assert(end_ - start_ == tpl_.size());
    assert(!pinStart_ || start_ == 0);
    // cannot test this unfortunately =(
    //   assert(!pinEnd_ || end_ == tpl_.size());
}

boost::optional<MutatedTemplate> Template::Mutate(const Mutation& mut)
{
    if (Length() == 0 && mut.LengthDiff() < 1) return boost::none;
    if (!InRange(mut.Start(), mut.End())) return boost::none;

    const uint8_t idx = detail::TranslationTable[static_cast<uint8_t>(mut.Base)];
    if (mut.Type != MutationType::DELETION && idx > 3)
        throw std::invalid_argument("invalid character in template!");

    return boost::make_optional(MutatedTemplate(*this, mut));
}

bool Template::ApplyMutation(const Mutation& mut)
{
    const uint8_t idx = detail::TranslationTable[static_cast<uint8_t>(mut.Base)];
    if (mut.Type != MutationType::DELETION && idx > 3)
        throw std::invalid_argument("invalid character in template!");

    bool mutApplied = false;

    if (Length() == 0 && mut.LengthDiff() < 1) goto finish;
    if (!InRange(mut.Start(), mut.End())) goto finish;

    {
        const size_t i = mut.Start() - start_;

        if (mut.Type == MutationType::INSERTION) {
            if (i > 0) tpl_[i - 1] = cfg_->Populate({tpl_[i - 1].Base, mut.Base})[0];

            if (i < tpl_.size())
                tpl_.insert(tpl_.begin() + i, cfg_->Populate({mut.Base, tpl_[i].Base})[0]);
            else
                tpl_.emplace_back(TemplatePosition{mut.Base, idx, 1.0, 0.0, 0.0, 0.0});
        } else if (mut.Type == MutationType::SUBSTITUTION) {
            if (i > 0) tpl_[i - 1] = cfg_->Populate({tpl_[i - 1].Base, mut.Base})[0];

            if (i + 1 < tpl_.size())
                tpl_[i] = cfg_->Populate({mut.Base, tpl_[i + 1].Base})[0];
            else
                tpl_[i] = TemplatePosition{mut.Base, idx, 1.0, 0.0, 0.0, 0.0};
        } else if (mut.Type == MutationType::DELETION) {
            tpl_.erase(tpl_.begin() + i);

            if (i > 0) {
                if (i < tpl_.size())
                    tpl_[i - 1] = cfg_->Populate({tpl_[i - 1].Base, tpl_[i].Base})[0];
                else
                    tpl_[i - 1] =
                        TemplatePosition{tpl_[i - 1].Base, tpl_[i - 1].Idx, 1.0, 0.0, 0.0, 0.0};
            }
        } else
            throw std::invalid_argument(
                "invalid mutation type! must be DELETION, INSERTION, or "
                "SUBSTITUTION");

        mutApplied = true;
    }

finish:
    // update the start_ and end_ mappings
    AbstractTemplate::ApplyMutation(mut);

    assert(tpl_.size() == end_ - start_);
    assert(Length() == 0 ||
           ((*this)[Length() - 1].Match == 1.0 && (*this)[Length() - 1].Branch == 0.0 &&
            (*this)[Length() - 1].Stick == 0.0 && (*this)[Length() - 1].Deletion == 0.0));

    assert(!pinStart_ || start_ == 0);

    if (Length() < 2) throw TemplateTooSmall();

    return mutApplied;
}

size_t Template::Length() const { return tpl_.size(); }

const TemplatePosition& Template::operator[](size_t i) const { return tpl_[i]; }

std::unique_ptr<AbstractRecursor> Template::CreateRecursor(const PacBio::Data::MappedRead& mr,
                                                           double scoreDiff) const
{
    return cfg_->CreateRecursor(mr, scoreDiff);
}
double Template::ExpectedLLForEmission(MoveType move, uint8_t prev, uint8_t curr,
                                       MomentType moment) const
{
    return cfg_->ExpectedLLForEmission(move, prev, curr, moment);
}

//
// MutatedTemplate Function Definitions
//

MutatedTemplate::MutatedTemplate(const Template& master, const Mutation& mut)
    : AbstractTemplate(master.start_, master.end_, master.pinStart_, master.pinEnd_)
    , master_(master)
    , mut_(mut)
{
    // Sanity check the input arguments and calculate the relative position of our mutation
    assert(!pinStart_ || start_ == 0);
    assert(!pinEnd_ || end_ - start_ == master_.tpl_.size());
    mutStart_ = mut.Start() - master_.start_;
    mutEnd_ = mut.End() - master_.start_;

    // Out-of-range idx values should be caught before construction in Template::Mutate()
    const uint8_t idx = detail::TranslationTable[static_cast<uint8_t>(mut.Base)];

    // Depending on the MutationType, fill out the 2-base MutTpl_ with the model parameters for
    // the base where the mutation occurs and (if possible) the base that succeeds it to account
    // for the context change caused by the mutation
    //
    // All mutations below are described as Before(B), Position(P), or After(A) with Mutated(M)
    // for the new nucleotide, such that the pre-mutation template can be read as "B-P-A".  Since
    // the parameters of all bases are induced by the context of the preceeding base, the context
    // of the targeted position is "B-P", and the context of the successor base is "P-A".
    if (mut.Type == MutationType::INSERTION) {
        // For Insertions the new sequence is "B-M-P-A", the mutation context is now "B-M" and
        // the context of the successor position is "M-P"
        if (mutStart_ > 0)
            mutTpl_[0] = master_.cfg_->Populate({master_.tpl_[mutStart_ - 1].Base, mut.Base})[0];

        if (mutStart_ < master_.tpl_.size())
            mutTpl_[1] = master_.cfg_->Populate({mut.Base, master_.tpl_[mutStart_].Base})[0];
        else
            mutTpl_[1] = TemplatePosition{mut.Base, idx, 1.0, 0.0, 0.0, 0.0};
    } else if (mut.Type == MutationType::SUBSTITUTION) {
        // For Substitutions the new sequence is "B-M-A", the mutation context is now "B-M" and
        // the context of the successor position is "M-A"
        if (mutStart_ > 0)
            mutTpl_[0] = master_.cfg_->Populate({master_.tpl_[mutStart_ - 1].Base, mut.Base})[0];

        if (mutStart_ + 1 < master_.tpl_.size())
            mutTpl_[1] = master_.cfg_->Populate({mut.Base, master_.tpl_[mutStart_ + 1].Base})[0];
        else
            mutTpl_[1] = TemplatePosition{mut.Base, idx, 1.0, 0.0, 0.0, 0.0};
    } else if (mut.Type == MutationType::DELETION) {
        // For Deletions the new sequence is "B-A", the mutational context is now "B-A" and
        // the context of the successor position is whatever "A"s context was previously
        if (mutStart_ > 0) {
            if (mutStart_ + 1 < master_.tpl_.size())
                mutTpl_[0] = master_.cfg_->Populate(
                    {master_.tpl_[mutStart_ - 1].Base, master_.tpl_[mutStart_ + 1].Base})[0];
            else
                mutTpl_[0] = TemplatePosition{master.tpl_[mutStart_ - 1].Base,
                                              master.tpl_[mutStart_ - 1].Idx,
                                              1.0,
                                              0.0,
                                              0.0,
                                              0.0};
        }

        // if there's a successor, the params for the mutant position
        //   are equiv to the next position in the original tpl
        if (mutStart_ + 1 < master_.tpl_.size()) mutTpl_[1] = master_.tpl_[mutStart_ + 1];
    } else
        throw std::invalid_argument(
            "invalid mutation type! must be DELETION, INSERTION, or "
            "SUBSTITUTION");

    mutOff_ = mut.LengthDiff();

    assert(Length() == 0 ||
           ((*this)[Length() - 1].Match == 1.0 && (*this)[Length() - 1].Branch == 0.0 &&
            (*this)[Length() - 1].Stick == 0.0 && (*this)[Length() - 1].Deletion == 0.0));
}

bool MutatedTemplate::ApplyMutation(const Mutation& mut)
{
    assert(!pinStart_ || start_ == 0);
    const bool mutApplied = AbstractTemplate::ApplyMutation(mut);
    assert(!pinStart_ || start_ == 0);
    return mutApplied;
}

size_t MutatedTemplate::Length() const { return end_ - start_ + mutOff_; }

size_t MutatedTemplate::MutationStart() const { return mutStart_; }

size_t MutatedTemplate::MutationEnd() const
{
    if (mut_.Type == MutationType::INSERTION || mut_.Type == MutationType::ANY_INSERTION)
        return mutStart_;

    // if (mut_.Type == MutationType::SUBSTITUTION ||
    //     mut_.Type == MutationType::ANY_SUBSTITUTION ||
    //     mut_.Type == MutationType::DELETION)
    return mutStart_ + 1;
}

int MutatedTemplate::LengthDiff() const { return mut_.LengthDiff(); }

const TemplatePosition& MutatedTemplate::operator[](const size_t i) const
{
    // For everything up to the base before mutStart_, just return what we have
    if (i + 1 < mutStart_) return master_.tpl_[i];

    // if we're beyond the mutation position, we have to adjust for any change in
    // template length caused by the mutation before returning
    else if (i > mutStart_)
        return master_.tpl_[i - mutOff_];

    // otherwise if we're the base before mutStart_, 0, else 1 of our mutated tpl
    // params
    return mutTpl_[i == mutStart_];
}

boost::optional<MutatedTemplate> MutatedTemplate::Mutate(const Mutation& mut)
{
    return boost::none;
}

std::unique_ptr<AbstractRecursor> MutatedTemplate::CreateRecursor(
    const PacBio::Data::MappedRead& mr, double scoreDiff) const
{
    return master_.CreateRecursor(mr, scoreDiff);
}

double MutatedTemplate::ExpectedLLForEmission(MoveType move, uint8_t prev, uint8_t curr,
                                              MomentType moment) const
{
    return master_.ExpectedLLForEmission(move, prev, curr, moment);
}

}  // namespace Consensus
}  // namespace PacBio
