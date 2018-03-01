// Author: Lance Hepler

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

#include <boost/optional.hpp>

#include <pacbio/align/LinearAlignment.h>
#include <pacbio/exception/InvalidEvaluatorException.h>

#include "Constants.h"
#include "EvaluatorImpl.h"
#include "matrix/BasicDenseMatrix.h"

using namespace PacBio::Data;
using namespace PacBio::Exception;

namespace PacBio {
namespace Consensus {
namespace {  // anonymous

static constexpr const double ALPHA_BETA_MISMATCH_TOLERANCE = 0.001;
static constexpr const double EARLY_ALPHA_BETA_MISMATCH_TOLERANCE = 0.0001;

#if 0
std::ostream& operator<<(std::ostream& out, const std::pair<size_t, size_t>& x)
{
    return out << '(' << x.first << ", " << x.second << ')';
}

void WriteMatrix(const ScaledMatrix& mat)
{
    std::cerr << std::pair<size_t, size_t>(mat.Rows(), mat.Columns()) << std::endl;

    for (size_t j = 0; j < mat.Columns(); ++j)
        std::cerr << " " << mat.UsedRowRange(j);
    std::cerr << std::endl;

    std::cerr << "lg: ";
    for (size_t j = 0; j < mat.Columns(); ++j)
        std::cerr << "\t" << std::fixed << std::setprecision(3) << mat.GetLogScale(j);
    std::cerr << std::endl;

    std::cerr << "lgS: ";
    double lgS = 0.0;
    for (size_t j = 0; j < mat.Columns(); ++j)
        std::cerr << "\t" << std::fixed << std::setprecision(3) << (lgS += mat.GetLogScale(j));
    std::cerr << std::endl;

    for (size_t i = 0; i < mat.Rows(); ++i)
    {
        for (size_t j = 0; j < mat.Columns(); ++j)
        {
            std::cerr << "\t" << std::fixed << std::setprecision(3) << std::log(mat.Get(i, j)) + mat.GetLogScale(j);
        }
        std::cerr << std::endl;
    }
}
#endif

}  // namespace anonymous

EvaluatorImpl::EvaluatorImpl(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr,
                             const double scoreDiff)
    : tpl_{std::move(tpl)}
    , recursor_{tpl_->CreateRecursor(mr, scoreDiff)}
    , alpha_(mr.Length() + 1, tpl_->Length() + 1, ScaledMatrix::FORWARD)
    , beta_(mr.Length() + 1, tpl_->Length() + 1, ScaledMatrix::REVERSE)
    , extendBuffer_(mr.Length() + 1, EXTEND_BUFFER_COLUMNS, ScaledMatrix::FORWARD)
{
    numFlipFlops_ =
        recursor_->FillAlphaBeta(*tpl_, alpha_, beta_, EARLY_ALPHA_BETA_MISMATCH_TOLERANCE);
}

std::string EvaluatorImpl::ReadName() const { return recursor_->read_.Name; }

double EvaluatorImpl::LL(const Mutation& mut)
{
    // if we've masked out the mutation then just return the ll as-is
    if (mask_.Contains(mut)) return LL();

    // Make a View of the template of what it would look like w/ Mutation
    boost::optional<MutatedTemplate> mutTpl = tpl_->Mutate(mut);

    // if the mutation didn't hit this read, just return the ll as-is
    if (!mutTpl) return LL();

    // Otherwise calculate and return the score, modulo the CounterWeight
    size_t betaLinkCol = 1 + mutTpl->MutationEnd();
    size_t absoluteLinkColumn = 1 + mutTpl->MutationEnd() + mutTpl->LengthDiff();

    double score;

    bool atBegin = mutTpl->MutationStart() < 3;
    bool atEnd = (mutTpl->MutationEnd() + 3) > beta_.Columns();

    if (!atBegin && !atEnd) {
        const size_t extendLength = 2;
        const size_t extendStartCol = mutTpl->MutationStart() - mut.IsDeletion();

        extendBuffer_.SetDirection(ScaledMatrix::FORWARD);
        recursor_->ExtendAlpha(*mutTpl, alpha_, extendStartCol, extendBuffer_, extendLength);
        score = recursor_->LinkAlphaBeta(*mutTpl, extendBuffer_, extendLength, beta_, betaLinkCol,
                                         absoluteLinkColumn) +
                alpha_.GetLogProdScales(0, extendStartCol);
    } else if (!atBegin && atEnd) {
        //
        // Extend alpha to end
        //
        size_t extendStartCol = mutTpl->MutationStart() - 1;
        assert(mutTpl->Length() + 1 > extendStartCol);
        size_t extendLength = mutTpl->Length() - extendStartCol + 1;

        extendBuffer_.SetDirection(ScaledMatrix::FORWARD);
        recursor_->ExtendAlpha(*mutTpl, alpha_, extendStartCol, extendBuffer_, extendLength);
        score = std::log(extendBuffer_(recursor_->read_.Length(), extendLength - 1)) +
                alpha_.GetLogProdScales(0, extendStartCol) +
                extendBuffer_.GetLogProdScales(0, extendLength);
    } else if (atBegin && !atEnd) {
        // If the mutation occurs at positions 0 - 2
        size_t extendLastCol = mutTpl->MutationEnd();
        // We duplicate this math inside the function
        size_t extendLength = 1 + mutTpl->MutationEnd() + mutTpl->LengthDiff();

        extendBuffer_.SetDirection(ScaledMatrix::REVERSE);
        recursor_->ExtendBeta(*mutTpl, beta_, extendLastCol, extendBuffer_, mutTpl->LengthDiff());
        score = std::log(extendBuffer_(0, 0)) +
                beta_.GetLogProdScales(extendLastCol + 1, beta_.Columns()) +
                extendBuffer_.GetLogProdScales(0, extendLength);
    } else {
        assert(atBegin && atEnd);
        /* This should basically never happen...
           and is a total disaster if it does.  The basic idea is that
           FillAlpha and FillBeta use the real "template" while we test
           mutations using "virtual" template positions and the Extend/Link
           methods.  Trying to call FillAlpha to calculate the likelihood of a
        virtual
           mutation is therefore going to fail, as it calculates using the
           "real" template.
        throw TooSmallTemplateException();
         */

        //
        // Just do the whole fill
        //
        ScaledMatrix alphaP(recursor_->read_.Length() + 1, mutTpl->Length() + 1,
                            ScaledMatrix::FORWARD);
        recursor_->FillAlpha(*mutTpl, ScaledMatrix::Null(), alphaP);
        score = std::log(alphaP(recursor_->read_.Length(), mutTpl->Length())) +
                alphaP.GetLogProdScales();
    }

    return score + recursor_->UndoCounterWeights(recursor_->read_.Length());
}

double EvaluatorImpl::LL() const
{
    return std::log(beta_(0, 0)) + beta_.GetLogProdScales() +
           recursor_->UndoCounterWeights(recursor_->read_.Length());
}

std::pair<double, double> EvaluatorImpl::NormalParameters() const
{
    return tpl_->NormalParameters();
}

double EvaluatorImpl::ZScore() const
{
    double mean, var;
    std::tie(mean, var) = NormalParameters();
    return (LL() - mean) / std::sqrt(var);
}

inline void EvaluatorImpl::Recalculate()
{
    size_t I = recursor_->read_.Length() + 1;
    size_t J = tpl_->Length() + 1;
    alpha_.Reset(I, J);
    beta_.Reset(I, J);
    extendBuffer_.Reset(I, EXTEND_BUFFER_COLUMNS);
    recursor_->FillAlphaBeta(*tpl_, alpha_, beta_, ALPHA_BETA_MISMATCH_TOLERANCE);
}

bool EvaluatorImpl::ApplyMutation(const Mutation& mut)
{
    if (tpl_->ApplyMutation(mut)) {
        Recalculate();
        mask_.Mutate({mut});
        return true;
    }
    return false;
}

bool EvaluatorImpl::ApplyMutations(std::vector<Mutation>* muts)
{
    if (tpl_->ApplyMutations(muts)) {
        Recalculate();
        mask_.Mutate(*muts);
        return true;
    }
    return false;
}

const AbstractMatrix& EvaluatorImpl::Alpha() const { return alpha_; }

const AbstractMatrix& EvaluatorImpl::Beta() const { return beta_; }

const AbstractMatrix* EvaluatorImpl::AlphaView(MatrixViewConvention c) const
{
    auto* m = new BasicDenseMatrix(alpha_.Rows(), alpha_.Columns());

    for (size_t i = 0; i < alpha_.Rows(); ++i) {
        for (size_t j = 0; j < alpha_.Columns(); ++j) {
            switch (c) {
                case MatrixViewConvention::AS_IS:
                    (*m)(i, j) = alpha_(i, j);
                    break;
                case MatrixViewConvention::LOGSPACE:
                    (*m)(i, j) = std::log(alpha_(i, j)) + alpha_.GetLogScale(j);
                    break;
                case MatrixViewConvention::LOGPROBABILITY:
                    (*m)(i, j) = std::log(alpha_(i, j)) + alpha_.GetLogScale(j) +
                                 recursor_->UndoCounterWeights(i);
                    break;
            }
        }
    }

    return m;
}

const AbstractMatrix* EvaluatorImpl::BetaView(MatrixViewConvention c) const
{
    auto* m = new BasicDenseMatrix(beta_.Rows(), beta_.Columns());

    for (size_t i = 0; i < beta_.Rows(); ++i) {
        for (size_t j = 0; j < beta_.Columns(); ++j) {
            switch (c) {
                case MatrixViewConvention::AS_IS:
                    (*m)(i, j) = beta_(i, j);
                    break;
                case MatrixViewConvention::LOGSPACE:
                    (*m)(i, j) = std::log(beta_(i, j)) + beta_.GetLogScale(j);
                    break;
                case MatrixViewConvention::LOGPROBABILITY:
                    (*m)(i, j) = std::log(beta_(i, j)) + beta_.GetLogScale(j) +
                                 recursor_->UndoCounterWeights(beta_.Rows() - 1 - i);
                    break;
            }
        }
    }

    return m;
}

void EvaluatorImpl::MaskIntervals(const size_t radius, const double maxErrRate)
{
    using namespace PacBio::Align;

    const auto alignAndJustify = [this](const StrandType strand) {
        auto aln = std::unique_ptr<PairwiseAlignment>(
            AlignLinear(std::string(*tpl_), recursor_->read_.Seq));
        if (strand == StrandType::FORWARD)
            aln->Justify(LRType::LEFT);
        else if (strand == StrandType::REVERSE)
            aln->Justify(LRType::RIGHT);
        else if (strand == StrandType::UNMAPPED)
            throw InvalidEvaluatorException("Unmapped read in interval masking");
        else
            throw std::runtime_error("Unknown StrandType");
        return aln;
    };

    const auto aln = alignAndJustify(recursor_->read_.Strand);

    std::vector<size_t> errsBySite;
    size_t nErr = 0;
    for (const char c : aln->Transcript()) {
        switch (c) {
            case 'I':
                nErr += 1;
                break;
            case 'D':
            case 'R':
                nErr += 1;
            // fallthrough
            case 'M':
                errsBySite.emplace_back(nErr);
                nErr = 0;
                break;
            default:
                throw std::runtime_error("unknown value in Transcript");
        }
    }
    // terminal insertions
    if (!errsBySite.empty()) errsBySite.back() += nErr;

    if (errsBySite.size() != tpl_->Length()) throw std::runtime_error("|errsBySite| != |tpl|");

    // filter windows with extreme mutations
    const size_t start = tpl_->Start();
    for (size_t i = 0; i < errsBySite.size(); ++i) {
        const size_t b = (radius >= i) ? 0 : i - radius;
        const size_t e = std::min(i + radius + 1, errsBySite.size());
        nErr = 0;
        for (size_t j = b; j < e; ++j)
            nErr += errsBySite[j];
        const double errRate = static_cast<double>(nErr) / (e - b);
        if (errRate >= maxErrRate) mask_.Insert({start + b, start + e});
    }
}

}  // namespace Consensus
}  // namespace PacBio
