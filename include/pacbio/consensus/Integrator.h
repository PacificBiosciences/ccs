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

#include <cstdint>
#include <functional>
#include <iostream>
#include <memory>
#include <set>

#include <pacbio/consensus/Evaluator.h>
#include <pacbio/consensus/Exceptions.h>
#include <pacbio/consensus/Mutation.h>

namespace PacBio {
namespace Consensus {

struct IntegratorConfig
{
    double MinZScore;
    double ScoreDiff;

    IntegratorConfig(double minZScore = -3.5, double scoreDiff = 12.5);
};

enum struct AddReadResult : uint8_t
{
    SUCCESS,
    ALPHA_BETA_MISMATCH,
    POOR_ZSCORE,
    OTHER,
    /*
     This enum is used in other places to
     both size an array and index into it.  So
     users can know the number of elements needed in
     an array that could count valus of this enum, this
     should always be the last value in the enum.
     */
    SIZE
};

inline std::ostream& operator<<(std::ostream& os, AddReadResult result)
{
    static const char* names[] = {"SUCCESS", "ALPHA/BETA MISMATCH", "POOR Z-SCORE", "OTHER"};
    os << names[static_cast<size_t>(result)];
    return os;
}

std::set<std::string> SupportedChemistries();

class AbstractIntegrator
{
public:
    virtual ~AbstractIntegrator();

    virtual size_t TemplateLength() const = 0;

    virtual char operator[](size_t i) const = 0;
    virtual operator std::string() const = 0;

    virtual double LL(const Mutation& mut);
    virtual double LL() const;

    double AvgZScore() const;
    std::vector<double> ZScores() const;
    std::vector<std::pair<double, double>> NormalParameters() const;

    virtual void ApplyMutation(const Mutation& mut) = 0;
    virtual void ApplyMutations(std::vector<Mutation>* muts) = 0;

    virtual AddReadResult AddRead(const MappedRead& read) = 0;

    // For debugging purposes
    // (Note that these include results include all evaluators, even the inactive ones)
    std::vector<double> LLs(const Mutation& mut);
    std::vector<double> LLs() const;
    std::vector<std::string> ReadNames() const;

protected:
    Mutation ReverseComplement(const Mutation& mut) const;

    AbstractIntegrator(const IntegratorConfig& cfg);

    // move constructor
    AbstractIntegrator(AbstractIntegrator&&);

    AddReadResult AddRead(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& read);

    IntegratorConfig cfg_;
    std::vector<Evaluator> evals_;
};

class MonoMolecularIntegrator : public AbstractIntegrator
{
public:
    MonoMolecularIntegrator(const std::string& tpl, const IntegratorConfig& cfg, const SNR& snr,
                            const std::string& model);

    // move constructor
    MonoMolecularIntegrator(MonoMolecularIntegrator&&);

    size_t TemplateLength() const;

    char operator[](size_t i) const;
    operator std::string() const;

    double LL(const Mutation& mut);
    inline double LL() const { return AbstractIntegrator::LL(); }
    void ApplyMutation(const Mutation& mut);
    void ApplyMutations(std::vector<Mutation>* muts);

    AddReadResult AddRead(const MappedRead& read);

protected:
    std::string mdl_;
    SNR snr_;
    Template fwdTpl_;
    Template revTpl_;
};

class MultiMolecularIntegrator : public AbstractIntegrator
{
public:
    MultiMolecularIntegrator(const std::string& tpl, const IntegratorConfig& cfg);

    size_t TemplateLength() const;

    char operator[](size_t i) const;
    operator std::string() const;

    void ApplyMutation(const Mutation& mut);
    void ApplyMutations(std::vector<Mutation>* muts);

    AddReadResult AddRead(const MappedRead& read);

protected:
    std::unique_ptr<AbstractTemplate> GetTemplate(const MappedRead& read, const SNR& snr);

    std::string fwdTpl_;
    std::string revTpl_;

private:
    friend struct std::hash<MultiMolecularIntegrator>;
};

}  // namespace Consensus
}  // namespace PacBio
