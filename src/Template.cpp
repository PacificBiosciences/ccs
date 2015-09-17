
#include <cassert>
#include <cmath>
#include <stdexcept>

#include <pacbio/consensus/Template.h>

namespace PacBio {
namespace Consensus {

AbstractTemplate::~AbstractTemplate() {}
void AbstractTemplate::ApplyMutations(std::vector<Mutation>* const muts)
{
    // make sure the mutations are sorted by site: End() then Start()
    std::sort(muts->begin(), muts->end(), Mutation::SiteComparer);

    for (auto it = muts->crbegin(); it != muts->crend(); ++it)
        ApplyMutation(*it);
}

std::tuple<double, double> AbstractTemplate::NormalParameters(const size_t start,
                                                              const size_t end) const
{
    double mean = 0.0, var = 0.0;

    for (size_t i = start; i < end; ++i) {
        double m, v;
        std::tie(m, v) = SiteNormalParameters(i);
        mean += m;
        var += v;
    }

    return std::make_tuple(mean, var);
}

std::tuple<double, double> AbstractTemplate::SiteNormalParameters(const size_t i) const
{
    // std::log(1.0/3);
    constexpr double lgThird = -1.0986122886681098;

    const auto params = (*this)[i];

    const double p_m = params.Match, l_m = std::log(p_m), l2_m = l_m * l_m;
    const double p_d = params.Deletion, l_d = std::log(p_d), l2_d = l_d * l_d;
    const double p_b = params.Branch, l_b = std::log(p_b), l2_b = l_b * l_b;
    const double p_s = params.Stick, l_s = std::log(p_s), l2_s = l_s * l_s;

    const double eps = 1.0 - BaseEmissionPr(MoveType::MATCH, params.Base, params.Base);
    const double E_M = (1.0 - eps) * 0.0 + eps * lgThird, E2_M = eps * lgThird * lgThird;
    const double E_D = 0.0, E2_D = E_D * E_D;
    const double E_B = 0.0, E2_B = E_B * E_B;
    const double E_S = lgThird, E2_S = E_S * E_S;

    auto ENN = [=](const double l_m, const double l_d, const double l_b, const double l_s,
                   const double E_M, const double E_D, const double E_B, const double E_S) {
        const double E_MD = (l_m + E_M) * p_m / (p_m + p_d) + (l_d + E_D) * p_d / (p_m + p_d);
        const double E_I  = (l_b + E_B) * p_b / (p_b + p_s) + (l_s + E_S) * p_s / (p_b + p_s);
        const double E_BS = E_I * (p_s + p_b) / (p_m + p_d);
        return E_MD + E_BS;
    };

    const double mean = ENN(l_m, l_d, l_b, l_s, E_M, E_D, E_B, E_S);
    const double var  = ENN(l2_m, l2_d, l2_b, l2_s, E2_M, E2_D, E2_B, E2_S) - mean * mean;

    return std::make_tuple(mean, var);
}

Template::Template(const std::string& tpl, std::unique_ptr<ModelConfig>&& cfg)
    : cfg_(std::forward<std::unique_ptr<ModelConfig>>(cfg))
    , tpl_{cfg_->Populate(tpl)}
    , mutated_{false}
    , mutOff_{0}
{
}

size_t Template::Length() const { return tpl_.size() + mutOff_; }
const TemplatePosition& Template::operator[](size_t i) const
{
    // if no mutation, or everything up to the base before mutPos_, just return
    // what we have
    if (!IsMutated() || i + 1 < mutPos_) return tpl_[i];

    // if we're beyond the mutation position, take the appropriate base
    else if (i > mutPos_)
        return tpl_[i - mutOff_];

    // otherwise if we're the base before mutPos_, 0, else 1 of our mutated tpl
    // params
    return mutTpl_[i == mutPos_];
}

bool Template::IsMutated() const { return mutated_; }
void Template::Mutate(const Mutation& mut)
{
    Reset();
    mutPos_ = mut.Start();
    if (mut.Type == MutationType::INSERTION) {
        if (mutPos_ > 0) mutTpl_[0] = cfg_->Populate({tpl_[mutPos_ - 1].Base, mut.Base})[0];

        if (mutPos_ < tpl_.size())
            mutTpl_[1] = cfg_->Populate({mut.Base, tpl_[mutPos_].Base})[0];
        else
            mutTpl_[1] = TemplatePosition{mut.Base, 0.0, 0.0, 0.0, 0.0};
    } else if (mut.Type == MutationType::SUBSTITUTION) {
        if (mutPos_ > 0) mutTpl_[0] = cfg_->Populate({tpl_[mutPos_ - 1].Base, mut.Base})[0];

        if (mutPos_ + 1 < tpl_.size())
            mutTpl_[1] = cfg_->Populate({mut.Base, tpl_[mutPos_ + 1].Base})[0];
        else
            mutTpl_[1] = TemplatePosition{mut.Base, 0.0, 0.0, 0.0, 0.0};
    } else if (mut.Type == MutationType::DELETION) {
        // if there's a predecessor, fill it
        if (mutPos_ > 0) {
            if (mutPos_ + 1 < tpl_.size())
                mutTpl_[0] = cfg_->Populate({tpl_[mutPos_ - 1].Base, tpl_[mutPos_ + 1].Base})[0];
            else
                mutTpl_[0] = TemplatePosition{tpl_[mutPos_ - 1].Base, 0.0, 0.0, 0.0, 0.0};
        }

        // if there's a successor, the params for the mutant position
        //   are equiv to the next position in the original tpl
        if (mutPos_ + 1 < tpl_.size()) mutTpl_[1] = tpl_[mutPos_ + 1];
    } else
        throw std::invalid_argument(
            "invalid mutation type! must be DELETION, INSERTION, or "
            "SUBSTITUTION");

    mutOff_  = mut.LengthDiff();
    mutated_ = true;

    assert((*this)[Length() - 1].Match == 0.0 && (*this)[Length() - 1].Branch == 0.0 &&
           (*this)[Length() - 1].Stick == 0.0 && (*this)[Length() - 1].Deletion == 0.0);
}

void Template::Reset()
{
    mutOff_  = 0;
    mutated_ = false;
}

void Template::ApplyMutation(const Mutation& mut)
{
    const size_t i = mut.Start();

    if (mut.Type == MutationType::INSERTION) {
        if (i > 0) tpl_[i - 1] = cfg_->Populate({tpl_[i - 1].Base, mut.Base})[0];

        if (i < tpl_.size())
            tpl_.insert(tpl_.begin() + i, cfg_->Populate({mut.Base, tpl_[i].Base})[0]);
        else
            tpl_.emplace_back(TemplatePosition{mut.Base, 0.0, 0.0, 0.0, 0.0});
    } else if (mut.Type == MutationType::SUBSTITUTION) {
        if (i > 0) tpl_[i - 1] = cfg_->Populate({tpl_[i - 1].Base, mut.Base})[0];

        if (i + 1 < tpl_.size())
            tpl_[i] = cfg_->Populate({mut.Base, tpl_[i + 1].Base})[0];
        else
            tpl_[i].Base = mut.Base;
    } else if (mut.Type == MutationType::DELETION) {
        tpl_.erase(tpl_.begin() + i);

        if (i > 0) {
            if (i < tpl_.size())
                tpl_[i - 1] = cfg_->Populate({tpl_[i - 1].Base, tpl_[i].Base})[0];
            else
                tpl_[i - 1] = TemplatePosition{tpl_[i - 1].Base, 0.0, 0.0, 0.0, 0.0};
        }
    } else
        throw std::invalid_argument(
            "invalid mutation type! must be DELETION, INSERTION, or "
            "SUBSTITUTION");

    assert((*this)[Length() - 1].Match == 0.0 && (*this)[Length() - 1].Branch == 0.0 &&
           (*this)[Length() - 1].Stick == 0.0 && (*this)[Length() - 1].Deletion == 0.0);
}

VirtualTemplate::VirtualTemplate(const Template& master, size_t start, size_t end)
    : master_{master}, start_{start}, end_{end}
{
    assert(start_ < end_);
}

void VirtualTemplate::ApplyMutation(const Mutation& mut)
{
    // if the end of the mutation is before the end of our mapping,
    //   update the mapping
    if (mut.End() <= end_) end_ += mut.LengthDiff();

    // if the end of the mutation is before the start of our mapping,
    //   update the mapping
    if (mut.End() < start_) start_ += mut.LengthDiff();
}

}  // namespace Consensus
}  // namespace PacBio
