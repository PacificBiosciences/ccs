
#include <cassert>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include <pacbio/consensus/Template.h>

namespace PacBio {
namespace Consensus {

AbstractTemplate::AbstractTemplate(const size_t start, const size_t end, const bool pinStart,
                                   const bool pinEnd)
    : start_{start}, end_{end}, pinStart_{pinStart}, pinEnd_{pinEnd}
{
    assert(start_ <= end_);
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

size_t AbstractTemplate::TrueLength() const { return end_ - start_; }
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
    // TODO (ndelaney): This probably needs to be matched to the specific model being used.
    // std::log(1.0/3);
    constexpr double lgThird = -1.0986122886681098;
    const uint8_t prev = (i == 0) ? 0 : (*this)[i - 1].Idx;  // default base : A
    const auto params = (*this)[i];

    // get the substitution rate here:
    const double eps = SubstitutionRate(prev, params.Idx);

    const double p_m = params.Match, l_m = std::log(p_m), l2_m = l_m * l_m;
    const double p_d = params.Deletion, l_d = std::log(p_d), l2_d = l_d * l_d;
    const double p_b = params.Branch, l_b = std::log(p_b), l2_b = l_b * l_b;
    const double p_s = params.Stick, l_s = std::log(p_s), l2_s = l_s * l_s;

    const double lgeps = std::log(eps);
    const double lg1minusEps = std::log(1.0 - eps);

    // First moment expectations (zero terms used for clarity)
    const double E_M = (1.0 - eps) * lg1minusEps + eps * (lgThird + lgeps);
    const double E_D = 0.0;
    const double E_B = 0.0;
    const double E_S = lgThird;

    // Calculate first moment
    const double E_MD = (l_m + E_M) * p_m / (p_m + p_d) + (l_d + E_D) * p_d / (p_m + p_d);
    const double E_I = (l_b + E_B) * p_b / (p_b + p_s) + (l_s + E_S) * p_s / (p_b + p_s);
    const double E_BS = E_I * (p_s + p_b) / (p_m + p_d);
    const double mean = E_MD + E_BS;

    // Calculate second momment
    // Key expansion used repeatedly here: (A + B)^2 = A^2 + 2AB + B^2
    const double E2_M = (1.0 - eps) * pow(lg1minusEps, 2.0) + eps * pow(lgThird + lgeps, 2.0);
    const double E2_MD =
        (l2_m + 2 * l_m * E_M + E2_M) * p_m / (p_m + p_d) + l2_d * p_d / (p_m + p_d);
    const double E2_S = pow(lgThird, 2.0);
    const double E2_I =
        l2_b * p_b / (p_b + p_s) + (l2_s + 2 * E_S * l_s + E2_S) * p_s / (p_b + p_s);
    const double E2_BS = E2_I * (p_s + p_b) / (p_m + p_d);
    const double moment2 = E2_BS + 2 * E_BS * E_MD + E2_MD;
    const double var = moment2 - mean * mean;
    return std::make_pair(mean, var);
}

std::ostream& operator<<(std::ostream& os, const AbstractTemplate& tpl)
{
    os << std::string(tpl);
    return os;
}

Template::Template(const std::string& tpl, std::unique_ptr<ModelConfig>&& cfg)
    : Template(tpl, std::forward<std::unique_ptr<ModelConfig>>(cfg), 0, tpl.length(), true, true)
{
}

Template::Template(const std::string& tpl, std::unique_ptr<ModelConfig>&& cfg, const size_t start,
                   const size_t end, const bool pinStart, const bool pinEnd)
    : AbstractTemplate(start, end, pinStart, pinEnd)
    , cfg_(std::forward<std::unique_ptr<ModelConfig>>(cfg))
    , tpl_{cfg_->Populate(tpl)}
    , mutated_{false}
    , mutOff_{0}
{
    assert(end_ - start_ == tpl_.size());
    assert(!pinStart_ || start_ == 0);
    // cannot test this unfortunately =(
    //   assert(!pinEnd_ || end_ == tpl_.size());
}

boost::optional<Mutation> Template::Mutate(const Mutation& mut)
{
    Reset();

    if (Length() == 0 && mut.LengthDiff() < 1) return boost::none;
    if (!InRange(mut.Start(), mut.End())) return boost::none;

    const uint8_t idx = detail::TranslationTable[static_cast<uint8_t>(mut.Base)];

    mutStart_ = mut.Start() - start_;
    mutEnd_ = mut.End() - start_;

    if (mut.Type == MutationType::INSERTION) {
        if (mutStart_ > 0) mutTpl_[0] = cfg_->Populate({tpl_[mutStart_ - 1].Base, mut.Base})[0];

        if (mutStart_ < tpl_.size())
            mutTpl_[1] = cfg_->Populate({mut.Base, tpl_[mutStart_].Base})[0];
        else
            mutTpl_[1] = TemplatePosition{mut.Base, idx, 1.0, 0.0, 0.0, 0.0};
    } else if (mut.Type == MutationType::SUBSTITUTION) {
        if (mutStart_ > 0) mutTpl_[0] = cfg_->Populate({tpl_[mutStart_ - 1].Base, mut.Base})[0];

        if (mutStart_ + 1 < tpl_.size())
            mutTpl_[1] = cfg_->Populate({mut.Base, tpl_[mutStart_ + 1].Base})[0];
        else
            mutTpl_[1] = TemplatePosition{mut.Base, idx, 1.0, 0.0, 0.0, 0.0};
    } else if (mut.Type == MutationType::DELETION) {
        // if there's a predecessor, fill it
        if (mutStart_ > 0) {
            if (mutStart_ + 1 < tpl_.size())
                mutTpl_[0] =
                    cfg_->Populate({tpl_[mutStart_ - 1].Base, tpl_[mutStart_ + 1].Base})[0];
            else
                mutTpl_[0] = TemplatePosition{
                    tpl_[mutStart_ - 1].Base, tpl_[mutStart_ - 1].Idx, 1.0, 0.0, 0.0, 0.0};
        }

        // if there's a successor, the params for the mutant position
        //   are equiv to the next position in the original tpl
        if (mutStart_ + 1 < tpl_.size()) mutTpl_[1] = tpl_[mutStart_ + 1];
    } else
        throw std::invalid_argument(
            "invalid mutation type! must be DELETION, INSERTION, or "
            "SUBSTITUTION");

    mutOff_ = mut.LengthDiff();
    mutated_ = true;

    assert(Length() == 0 ||
           ((*this)[Length() - 1].Match == 1.0 && (*this)[Length() - 1].Branch == 0.0 &&
            (*this)[Length() - 1].Stick == 0.0 && (*this)[Length() - 1].Deletion == 0.0));

    return Mutation(mut.Type, mutStart_, mut.Base);
}

void Template::Reset()
{
    mutOff_ = 0;
    mutated_ = false;
}

bool Template::ApplyMutation(const Mutation& mut)
{
    if (mutated_)
        throw std::runtime_error(
            "cannot ApplyMutation to an already mutated Template, call Reset first");

    const uint8_t idx = detail::TranslationTable[static_cast<uint8_t>(mut.Base)];
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

    return mutApplied;
}

VirtualTemplate::VirtualTemplate(const Template& master, const size_t start, const size_t end,
                                 const bool pinStart, const bool pinEnd)
    : AbstractTemplate(start, end, pinStart, pinEnd), master_(master)
{
    assert(!pinStart_ || start_ == 0);
    assert(!pinEnd_ || end_ == master_.tpl_.size());
}

bool VirtualTemplate::ApplyMutation(const Mutation& mut)
{
    assert(!pinStart_ || start_ == 0);
    const bool mutApplied = AbstractTemplate::ApplyMutation(mut);
    assert(!pinStart_ || start_ == 0);
    return mutApplied;
}

}  // namespace Consensus
}  // namespace PacBio
