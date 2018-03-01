// Author: Lance Hepler

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

AbstractTemplate::~AbstractTemplate() = default;
AbstractTemplate::operator std::string() const
{
    std::ostringstream ss;
    for (size_t i = 0; i < Length(); ++i)
        ss << (*this)[i].Base;
    return ss.str();
}

boost::optional<MutatedTemplate> AbstractTemplate::Mutate(const Mutation& mut)
{
    if (Length() == 0 && mut.LengthDiff() < 1) return boost::none;
    if (!InRange(mut.Start(), mut.End())) return boost::none;
    const boost::optional<Mutation> tMut = mut.Translate(start_, Length());
    if (!tMut) return boost::none;
    return boost::make_optional(MutatedTemplate(*this, *tMut));
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
    const auto prev =
        (i == 0) ? AlleleRep::FromASCII('A') : (*this)[i - 1].Idx;  // default base : A
    const auto curr = params.Idx;

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
    : Template(tpl, std::move(cfg), 0, tpl.length(), true, true)
{
}

Template::Template(const std::string& tpl, std::unique_ptr<ModelConfig>&& cfg, const size_t start,
                   const size_t end, const bool pinStart, const bool pinEnd)
    : AbstractTemplate(start, end, pinStart, pinEnd)
    , cfg_(std::move(cfg))
    , tpl_{cfg_->Populate(tpl)}
{
    assert(end_ - start_ == tpl_.size());
    assert(!pinStart_ || start_ == 0);
    // cannot test this unfortunately =(
    //   assert(!pinEnd_ || end_ == tpl_.size());
}

bool Template::ApplyMutation(const Mutation& mut)
{
    bool mutApplied = false;

    if (Length() == 0 && mut.LengthDiff() < 1) goto finish;
    if (!InRange(mut.Start(), mut.End())) goto finish;

    {
        const size_t b = mut.Start() - start_;

        if (mut.IsDeletion()) {
            const size_t e = mut.End() - start_;
            tpl_.erase(tpl_.begin() + b, tpl_.begin() + e);

            if (b > 0) {
                if (b < tpl_.size())
                    tpl_[b - 1] = cfg_->Populate({tpl_[b - 1].Base, tpl_[b].Base})[0];
                else
                    tpl_[b - 1] = TemplatePosition{tpl_[b - 1].Base, 1.0, 0.0, 0.0, 0.0};
            }
        } else if (mut.IsInsertion()) {
            const auto elems = cfg_->Populate(mut.Bases());
            const size_t e = b + elems.size();

            tpl_.insert(tpl_.begin() + b, elems.begin(), elems.end());

            if (b > 0) tpl_[b - 1] = cfg_->Populate({tpl_[b - 1].Base, tpl_[b].Base})[0];
            if (0 < e && e < tpl_.size())
                tpl_[e - 1] = cfg_->Populate({tpl_[e - 1].Base, tpl_[e].Base})[0];
        } else if (mut.IsSubstitution()) {
            const auto elems = cfg_->Populate(mut.Bases());
            const size_t e = mut.End() - start_;

            for (size_t i = b; i < e; ++i)
                tpl_[i] = elems[i - b];

            if (b > 0) tpl_[b - 1] = cfg_->Populate({tpl_[b - 1].Base, tpl_[b].Base})[0];
            if (0 < e && e < tpl_.size())
                tpl_[e - 1] = cfg_->Populate({tpl_[e - 1].Base, tpl_[e].Base})[0];
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
double Template::ExpectedLLForEmission(MoveType move, const AlleleRep& prev, const AlleleRep& curr,
                                       MomentType moment) const
{
    return cfg_->ExpectedLLForEmission(move, prev, curr, moment);
}

//
// MutatedTemplate Function Definitions
//

// TODO(lhepler): this is most certainly busted for a deletion at 0 w/o pinStart_
MutatedTemplate::MutatedTemplate(const AbstractTemplate& master, const Mutation& mut)
    : AbstractTemplate(master.start_, master.end_, master.pinStart_, master.pinEnd_)
    , master_{master}
    , mut_{mut}
    , mutStart_{(mut.Start() > 0) ? mut.Start() - 1 : 0}
    , mutOff_{mut.LengthDiff()}
{
    // Sanity check the input arguments and calculate the relative position of our mutation
    assert(!pinStart_ || start_ == 0);
    assert(!pinEnd_ || end_ - start_ == master_.Length());

    // Fill out mutTpl_ with the model parameters for the base before the mutation and all the
    // bases changed by it.
    //
    // All mutations below are described as Before(B), Position(P), or After(A) with Mutated(M)
    // for the new nucleotide, such that the pre-mutation template can be read as "B-P-A".  Since
    // the parameters of all bases are induced by the context of the preceeding base, the context
    // of the targeted position is "B-P", and the context of the successor base is "P-A".
    const size_t mStart = mut_.Start();
    const size_t mEnd = mut_.End();
    if (mut_.IsDeletion()) {
        if (mStart > 0) {
            if (mEnd < master_.Length()) {
                mutTpl_.emplace_back(
                    master_.Config().Populate({master_[mStart - 1].Base, master_[mEnd].Base})[0]);
            } else
                mutTpl_.emplace_back(TemplatePosition{master[mStart - 1].Base, 1.0, 0.0, 0.0, 0.0});
        }
    } else if (mut_.IsInsertion() || mut_.IsSubstitution()) {
        if (mStart > 0)
            mutTpl_.emplace_back(
                master_.Config().Populate({master_[mStart - 1].Base, mut_.Bases()[0]})[0]);

        auto elems = master_.Config().Populate(mut.Bases());
        mutTpl_.reserve(mutTpl_.size() + elems.size());
        mutTpl_.insert(mutTpl_.end(), elems.begin(), elems.end());

        if (mEnd < master_.Length())
            mutTpl_.back() =
                master_.Config().Populate({mut_.Bases().back(), master_[mEnd].Base})[0];
    } else
        throw std::invalid_argument(
            "invalid mutation type! must be DELETION, INSERTION, or "
            "SUBSTITUTION");

    assert((mut_.IsDeletion() && mutTpl_.size() == (mut.Start() > 0)) ||
           (mutTpl_.size() == mut_.Bases().size() + (mut.Start() > 0)));
    assert(Length() == 0 ||
           ((*this)[Length() - 1].Match == 1.0 && (*this)[Length() - 1].Branch == 0.0 &&
            (*this)[Length() - 1].Stick == 0.0 && (*this)[Length() - 1].Deletion == 0.0));
}

bool MutatedTemplate::ApplyMutation(const Mutation& mut)
{
    throw std::runtime_error("MutatedTemplate cannot perform ApplyMutation!");
}

size_t MutatedTemplate::Length() const { return end_ - start_ + mutOff_; }

MutationType MutatedTemplate::Type() const { return mut_.Type(); }

size_t MutatedTemplate::MutationStart() const { return mut_.Start(); }

size_t MutatedTemplate::MutationEnd() const { return mut_.End(); }

int MutatedTemplate::LengthDiff() const { return mutOff_; }

const TemplatePosition& MutatedTemplate::operator[](const size_t i) const
{
    // For everything up to the base before mutStart_, just return what we have
    if (i < mutStart_) return master_[i];

    // if we're beyond the mutation position, we have to adjust for any change in
    // template length caused by the mutation before returning
    else if (i >= mutStart_ + mutTpl_.size())
        return master_[i - mutOff_];

    return mutTpl_[i - mutStart_];
}

std::unique_ptr<AbstractRecursor> MutatedTemplate::CreateRecursor(
    const PacBio::Data::MappedRead& mr, double scoreDiff) const
{
    return master_.CreateRecursor(mr, scoreDiff);
}

double MutatedTemplate::ExpectedLLForEmission(MoveType move, const AlleleRep& prev,
                                              const AlleleRep& curr, MomentType moment) const
{
    return master_.ExpectedLLForEmission(move, prev, curr, moment);
}

}  // namespace Consensus
}  // namespace PacBio
