// Copyright (c) 2011-2014, Pacific Biosciences of California, Inc.
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

// Author: David Alexander

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <iostream>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <pacbio/consensus/ModelConfig.h>
#include <pacbio/consensus/Mutation.h>
#include <pacbio/consensus/Read.h>
#include <pacbio/consensus/Template.h>

#include "../src/ModelFactory.h"
#include "Mutations.h"
#include "RandomDNA.h"

using namespace std;
using namespace PacBio::Consensus;  // NOLINT

using ::testing::UnorderedElementsAreArray;

namespace PacBio {
namespace Consensus {

bool operator==(const TemplatePosition& lhs, const TemplatePosition& rhs)
{
    return lhs.Base == rhs.Base && lhs.Match == rhs.Match && lhs.Branch == rhs.Branch &&
           lhs.Stick == rhs.Stick && lhs.Deletion == rhs.Deletion;
}

bool operator!=(const TemplatePosition& lhs, const TemplatePosition& rhs) { return !(lhs == rhs); }
bool operator==(const AbstractTemplate& lhs, const AbstractTemplate& rhs)
{
    if (lhs.Length() != rhs.Length()) return false;

    const size_t len = lhs.Length();

    // do not test the last base, can vary
    for (size_t i = 0; i + 1 < len; ++i)
        if (lhs[i] != rhs[i]) return false;

    // check the last base
    if (len > 0 && lhs[len - 1].Base != rhs[len - 1].Base) return false;

    return true;
}

bool operator!=(const AbstractTemplate& lhs, const AbstractTemplate& rhs) { return !(lhs == rhs); }
#if 0
ostream& operator<<(ostream& os, const AbstractTemplate& tpl)
{
    os << "Template:" << endl;
    for (size_t i = 0; i < tpl.Length(); ++i)
        os << "  " << tpl[i] << endl;
    return os;
}
#endif

}  // namespace Consensus
}  // namespace PacBio

namespace {

const string mdl = "P6-C4";
const SNR snr(10, 7, 5, 11);

string ApplyMutations(const string& tpl, vector<Mutation>&& muts)
{
    return ApplyMutations(tpl, &muts);
}

TEST(TemplateTest, ApplyMutations)
{
    constexpr auto INS = MutationType::INSERTION;
    constexpr auto DEL = MutationType::DELETION;
    constexpr auto SUB = MutationType::SUBSTITUTION;

    // insertion
    EXPECT_EQ("ACGT", ApplyMutations("CGT", {Mutation(INS, 0, 'A')}));
    EXPECT_EQ("ACGT", ApplyMutations("AGT", {Mutation(INS, 1, 'C')}));
    EXPECT_EQ("ACGT", ApplyMutations("ACT", {Mutation(INS, 2, 'G')}));
    EXPECT_EQ("ACGT", ApplyMutations("ACG", {Mutation(INS, 3, 'T')}));

    // substitution
    EXPECT_EQ("ACGT", ApplyMutations("XCGT", {Mutation(SUB, 0, 'A')}));
    EXPECT_EQ("ACGT", ApplyMutations("AXGT", {Mutation(SUB, 1, 'C')}));
    EXPECT_EQ("ACGT", ApplyMutations("ACXT", {Mutation(SUB, 2, 'G')}));
    EXPECT_EQ("ACGT", ApplyMutations("ACGX", {Mutation(SUB, 3, 'T')}));

    // deletion
    EXPECT_EQ("ACGT", ApplyMutations("XACGT", {Mutation(DEL, 0)}));
    EXPECT_EQ("ACGT", ApplyMutations("AXCGT", {Mutation(DEL, 1)}));
    EXPECT_EQ("ACGT", ApplyMutations("ACXGT", {Mutation(DEL, 2)}));
    EXPECT_EQ("ACGT", ApplyMutations("ACGXT", {Mutation(DEL, 3)}));
    EXPECT_EQ("ACGT", ApplyMutations("ACGTX", {Mutation(DEL, 4)}));
}

string ToString(const AbstractTemplate& tpl)
{
    string result;
    result.reserve(tpl.Length());
    for (size_t i = 0; i < tpl.Length(); ++i)
        result.push_back(tpl[i].Base);
    return result;
}

void TemplateEquivalence(const size_t nsamp, const size_t nvirt, const size_t len = 100)
{
    std::random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<size_t> randIdx(0, len - 1);
    bernoulli_distribution randPin(0.5);
    bernoulli_distribution randSpanning(0.33);

    for (size_t i = 0; i < nsamp; ++i) {
        const string tpl = RandomDNA(len, &gen);
        Template master(tpl, ModelFactory::Create(mdl, snr));
        EXPECT_EQ(tpl, ToString(master));
        vector<tuple<size_t, size_t, bool, bool>> coords;
        vector<VirtualTemplate> vtpls;
        vector<Template> rtpls;
        for (size_t j = 0; j < nvirt; ++j) {
            size_t start = 0;
            size_t end = len;
            // roughly 33% of the time having a spanning read
            if (randSpanning(gen)) {
                do {
                    start = randIdx(gen);
                    end = randIdx(gen);
                } while (start == end);
                if (end < start) swap(start, end);
                ++end;  // increment b by one (end-exclusive)
            }
            const bool pinStart = (start == 0) ? randPin(gen) : false;
            const bool pinEnd = (end == len) ? randPin(gen) : false;
            ASSERT_LT(start, end);
            coords.emplace_back(make_tuple(start, end, pinStart, pinEnd));
            vtpls.emplace_back(VirtualTemplate(master, start, end, pinStart, pinEnd));
            const string vtpl = tpl.substr(start, end - start);
            rtpls.emplace_back(
                Template(vtpl, ModelFactory::Create(mdl, snr), start, end, pinStart, pinEnd));
            EXPECT_EQ(vtpl, ToString(vtpls.back()));
            EXPECT_EQ(end - start, vtpls.back().Length());
        }
        vector<Mutation> mutations = Mutations(tpl);
        for (const auto& mut : mutations) {
            vector<Mutation> muts{mut};
            const string app = ApplyMutations(tpl, &muts);
            master.Mutate(mut);
            EXPECT_EQ(app, ToString(master));
            {
                Template mutated(app, ModelFactory::Create(mdl, snr));
                EXPECT_EQ(mutated, master);
            }
            for (size_t j = 0; j < nvirt; ++j) {
                size_t start, end;
                bool pinStart, pinEnd;
                tie(start, end, pinStart, pinEnd) = coords[j];
                const string vtpl = tpl.substr(start, end - start);
                vector<Mutation> vmuts;
                if ((pinStart || start < mut.End()) && (pinEnd || mut.Start() < end))
                    vmuts.emplace_back(Mutation(mut.Type, mut.Start() - start, mut.Base));
                const string vapp = ApplyMutations(vtpl, &vmuts);
                vtpls[j].Mutate(mut);
                rtpls[j].Mutate(mut);
                EXPECT_EQ(vapp, ToString(rtpls[j]));
                EXPECT_EQ(vapp, ToString(vtpls[j]));
                if (vapp != ToString(vtpls[j])) {
                    const size_t c =
                        start + ((pinStart || mut.End() < start) ? mut.LengthDiff() : 0);
                    const size_t d = end + ((pinEnd || mut.Start() < end) ? mut.LengthDiff() : 0);
                    cerr << "mut:  " << mut << endl;
                    if (!vmuts.empty()) cerr << "vmut: " << vmuts.back() << endl;
                    cerr << "off:  " << mut.LengthDiff() << endl;
                    cerr << "tpl:  " << tpl << endl;
                    cerr << "s,e:  " << start << "," << end << endl;
                    cerr << "ps,e: " << pinStart << "," << pinEnd << endl;
                    cerr << "ms,e: " << mut.Start() << "," << mut.End() << endl;
                    cerr << "app:  " << app << endl;
                    cerr << "c,d:  " << c << "," << d << endl;
                    cerr << "vtpl: " << vtpl << endl;
                    cerr << "vapp: " << vapp << endl;
                    // cerr << "sub:  " << app.substr(a + ((mut.Start() <= a) ? mut.LengthDiff() :
                    // 0), vapp.length()) << endl;
                    cerr << "vchl: " << ToString(vtpls[j]) << endl;
                    cerr << "rchl: " << ToString(rtpls[j]) << endl;
                    cerr << "vxxx: " << ToString(master).substr(c, d - c) << endl;
                }
                const int off = vmuts.empty() ? 0 : vmuts.back().LengthDiff();
                EXPECT_EQ(end - start + off, vtpls[j].Length());
                EXPECT_EQ(end - start + off, rtpls[j].Length());
                {
                    Template child(vapp, ModelFactory::Create(mdl, snr));
                    EXPECT_EQ(child, vtpls[j]);
                    EXPECT_EQ(child, rtpls[j]);
                }
                vtpls[j].Reset();
                rtpls[j].Reset();
            }
            master.Reset();
        }
    }
}

TEST(TemplateTest, TestVirtualTemplateEquivalence)
{
    TemplateEquivalence(1000, 20, 10);
    TemplateEquivalence(500, 20, 30);
}

TEST(TemplateTest, TestPinning)
{
    constexpr size_t len = 5;
    const string tpl(len, 'C');
    const string A(1, 'A');

    // pinStart and pinEnd
    {
        Template master(tpl, ModelFactory::Create(mdl, snr), 0, len, true, true);
        EXPECT_EQ(len, master.Length());
        EXPECT_EQ(tpl, ToString(master));
        master.ApplyMutation(Mutation(MutationType::INSERTION, len, 'A'));
        EXPECT_EQ(len + 1, master.Length());
        EXPECT_EQ(tpl + A, ToString(master));
        master.ApplyMutation(Mutation(MutationType::INSERTION, 0, 'A'));
        EXPECT_EQ(len + 2, master.Length());
        EXPECT_EQ(A + tpl + A, ToString(master));
    }
    // no pinStart but pinEnd
    {
        Template master(tpl, ModelFactory::Create(mdl, snr), 0, len, false, true);
        EXPECT_EQ(len, master.Length());
        EXPECT_EQ(tpl, ToString(master));
        master.ApplyMutation(Mutation(MutationType::INSERTION, 0, 'A'));
        EXPECT_EQ(len, master.Length());
        EXPECT_EQ(tpl, ToString(master));
        // the coords are now 1..6, so a new terminal mutation is at len+1
        master.ApplyMutation(Mutation(MutationType::INSERTION, len + 1, 'A'));
        EXPECT_EQ(len + 1, master.Length());
        EXPECT_EQ(tpl + A, ToString(master));
    }
    // pinStart but no pinEnd
    {
        Template master(tpl, ModelFactory::Create(mdl, snr), 0, len, true, false);
        EXPECT_EQ(len, master.Length());
        EXPECT_EQ(tpl, ToString(master));
        master.ApplyMutation(Mutation(MutationType::INSERTION, len, 'A'));
        EXPECT_EQ(len, master.Length());
        EXPECT_EQ(tpl, ToString(master));
        master.ApplyMutation(Mutation(MutationType::INSERTION, 0, 'A'));
        EXPECT_EQ(len + 1, master.Length());
        EXPECT_EQ(A + tpl, ToString(master));
    }
    // no pinStart or pinEnd
    {
        Template master(tpl, ModelFactory::Create(mdl, snr), 0, len, false, false);
        EXPECT_EQ(len, master.Length());
        EXPECT_EQ(tpl, ToString(master));
        master.ApplyMutation(Mutation(MutationType::INSERTION, len, 'A'));
        EXPECT_EQ(len, master.Length());
        EXPECT_EQ(tpl, ToString(master));
        master.ApplyMutation(Mutation(MutationType::INSERTION, 0, 'A'));
        EXPECT_EQ(len, master.Length());
        EXPECT_EQ(tpl, ToString(master));
    }
}

TEST(TemplateTest, NullTemplate)
{
    const string tpl = "ACGT";
    const size_t len = tpl.length();
    const Mutation mut(MutationType::DELETION, 0, '-');

    Template master(tpl, ModelFactory::Create(mdl, snr), 0, len, true, true);
    VirtualTemplate virt(master, 0, len, false, false);

    EXPECT_EQ(len, master.Length());

    for (size_t i = 1; i <= len; ++i) {
        master.ApplyMutation(mut);
        EXPECT_EQ(len - i, master.Length());
        virt.ApplyMutation(mut);
        EXPECT_EQ(len - i, virt.Length());
    }

    {
        const string A(1, 'A');

        auto mut_ = master.Mutate(mut);
        EXPECT_FALSE(bool(mut_));
        master.ApplyMutation(mut);

        mut_ = master.Mutate(Mutation(MutationType::INSERTION, 0, 'A'));
        ASSERT_TRUE(bool(mut_));
        EXPECT_EQ(A, ToString(master));
        master.Reset();
        master.ApplyMutation(*mut_);
        EXPECT_EQ(A, ToString(master));
        virt.ApplyMutation(*mut_);
        EXPECT_EQ(0, virt.Length());
    }
}

}  // namespace anonymous
