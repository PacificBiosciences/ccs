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
#include <pacbio/consensus/Template.h>
#include <pacbio/data/Read.h>
#include <pacbio/exception/StateError.h>

#include "../src/ModelFactory.h"
#include "Mutations.h"
#include "RandomDNA.h"

using std::string;
using std::tuple;
using std::vector;

using std::cerr;
using std::endl;

using std::make_tuple;
using std::swap;
using std::tie;

using namespace PacBio::Consensus;  // NOLINT
using namespace PacBio::Exception;  // NOLINT

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

    // do not test the last base, can vary
    for (size_t i = 0; i < lhs.Length(); ++i)
        if (lhs[i] != rhs[i]) return false;

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

namespace TemplateTests {

const string mdl = "P6-C4";
const SNR snr(10, 7, 5, 11);

string ApplyMutations(const string& tpl, vector<Mutation>&& muts)
{
    return ApplyMutations(tpl, &muts);
}

TEST(TemplateTest, ApplyMutations)
{
    // insertion
    EXPECT_EQ("ACGT", ApplyMutations("CGT", {Mutation::Insertion(0, 'A')}));
    EXPECT_EQ("ACGT", ApplyMutations("AGT", {Mutation::Insertion(1, 'C')}));
    EXPECT_EQ("ACGT", ApplyMutations("ACT", {Mutation::Insertion(2, 'G')}));
    EXPECT_EQ("ACGT", ApplyMutations("ACG", {Mutation::Insertion(3, 'T')}));

    // substitution
    EXPECT_EQ("ACGT", ApplyMutations("XCGT", {Mutation::Substitution(0, 'A')}));
    EXPECT_EQ("ACGT", ApplyMutations("AXGT", {Mutation::Substitution(1, 'C')}));
    EXPECT_EQ("ACGT", ApplyMutations("ACXT", {Mutation::Substitution(2, 'G')}));
    EXPECT_EQ("ACGT", ApplyMutations("ACGX", {Mutation::Substitution(3, 'T')}));

    // deletion
    EXPECT_EQ("ACGT", ApplyMutations("XACGT", {Mutation::Deletion(0, 1)}));
    EXPECT_EQ("ACGT", ApplyMutations("AXCGT", {Mutation::Deletion(1, 1)}));
    EXPECT_EQ("ACGT", ApplyMutations("ACXGT", {Mutation::Deletion(2, 1)}));
    EXPECT_EQ("ACGT", ApplyMutations("ACGXT", {Mutation::Deletion(3, 1)}));
    EXPECT_EQ("ACGT", ApplyMutations("ACGTX", {Mutation::Deletion(4, 1)}));
}

string ToString(const AbstractTemplate& tpl)
{
    string result;
    result.reserve(tpl.Length());
    for (size_t i = 0; i < tpl.Length(); ++i)
        result.push_back(tpl[i].Base);
    return result;
}

void TemplateEquivalence(const size_t nSamples, const size_t nReads, const size_t len = 100)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> randIdx(0, len - 1);
    std::bernoulli_distribution randPin(0.5);
    std::bernoulli_distribution randNonSpanning(0.33);

    for (size_t i = 0; i < nSamples; ++i) {
        const string tpl = RandomDNA(len, &gen);               // "Reference" sequence
        Template master(tpl, ModelFactory::Create(mdl, snr));  // "Reference" template object
        EXPECT_EQ(tpl, ToString(master));

        // Generate a random mixture of spanning and non-spanning reads,
        //  as well as pinned and non-pinned reads for testing
        vector<tuple<size_t, size_t, bool, bool>> coords;  // "Read" coordinates
        vector<Template> rtpls;                            // "Read" templates
        for (size_t j = 0; j < nReads; ++j) {
            int start = 0;
            int end = len;
            // Roughly 33% of the time having a non-spanning read and pick
            //  start and end positions at random (minimum 2bp)
            if (randNonSpanning(gen)) {
                do {
                    start = randIdx(gen);
                    end = randIdx(gen);
                } while (std::abs(start - end) < 2);
                if (end < start) swap(start, end);
                ++end;  // increment b by one (end-exclusive)
            }
            ASSERT_LT(start, end);

            // Save both the resulting "read" template and it's positions
            const bool pinStart = (start == 0) ? randPin(gen) : false;
            const bool pinEnd = (end == static_cast<int>(len)) ? randPin(gen) : false;
            coords.emplace_back(make_tuple(start, end, pinStart, pinEnd));
            const string rtpl = tpl.substr(start, end - start);
            rtpls.emplace_back(
                Template(rtpl, ModelFactory::Create(mdl, snr), start, end, pinStart, pinEnd));
            EXPECT_EQ(end - start, rtpls.back().Length());
        }

        vector<Mutation> mutations = Mutations(tpl);
        for (const auto& mut : mutations) {

            // Applying a mutation to a template with Mutate() should produce
            //  the same string as the operation applied to the underlying string
            vector<Mutation> muts{mut};
            const string refMutStr = ApplyMutations(tpl, &muts);  // mutated "Reference" string
            const auto refMutTpl = master.Mutate(mut);            // mutated "Reference" template
            EXPECT_TRUE(bool(refMutTpl));
            EXPECT_EQ(refMutStr, ToString(*refMutTpl));
            {
                const Template mutated(refMutStr, ModelFactory::Create(mdl, snr));
                Template applied(tpl, ModelFactory::Create(mdl, snr));
                applied.ApplyMutations(&muts);
                EXPECT_EQ(mutated, *refMutTpl);
                EXPECT_EQ(mutated, applied);
            }

            // Applying a mutation to a "subread" Template with Mutate() should produce
            //  the same string as the operation applied to the underlying string OR
            //  boost::none if it's out of range
            for (size_t j = 0; j < nReads; ++j) {
                size_t start, end;
                bool pinStart, pinEnd;
                tie(start, end, pinStart, pinEnd) = coords[j];
                const string rStr = tpl.substr(start, end - start);  // "Read" string

                const auto mTpl = rtpls[j].Mutate(mut);  // Mutated "Read" template
                if (!mTpl) {  // If Mutate() returned None, assert that the mutation is out of range
                    const bool isInRange =
                        ((pinStart || start < mut.End()) && (pinEnd || mut.Start() < end));

                    // Print a report if we fail to mutate the template, but are within range
                    if (isInRange) {
                        const size_t c =
                            start + ((pinStart || mut.End() < start) ? mut.LengthDiff() : 0);
                        const size_t d =
                            end + ((pinEnd || mut.Start() < end) ? mut.LengthDiff() : 0);
                        cerr << "Mut:          " << mut << endl;
                        cerr << "Off:          " << mut.LengthDiff() << endl;
                        cerr << "Ref:          " << tpl << endl;
                        cerr << "Start,End:    " << start << "," << end << endl;
                        cerr << "PinStart,End: " << pinStart << "," << pinEnd << endl;
                        cerr << "MutStart,End: " << mut.Start() << "," << mut.End() << endl;
                        cerr << "refMutStr:    " << refMutStr << endl;
                        cerr << "c,d:          " << c << "," << d << endl;
                        cerr << "refStr:       " << rStr << endl;
                        cerr << "refTpl:       " << ToString(master).substr(c, d - c) << endl;
                    }

                    EXPECT_FALSE(isInRange);

                } else {  // Otherwise it should be in-range of the "Read" template
                    const auto rMut = mut.Translate(start, end - start);
                    EXPECT_TRUE(bool(rMut));
                    vector<Mutation> rMuts{*rMut};                     // "Read"-space mutations
                    const string mStr = ApplyMutations(rStr, &rMuts);  // Mutated "Read" string

                    // Print a report if the mutated template isn't correct
                    if (mStr != ToString(*mTpl)) {
                        const size_t c =
                            start + ((pinStart || mut.End() < start) ? mut.LengthDiff() : 0);
                        const size_t d =
                            end + ((pinEnd || mut.Start() < end) ? mut.LengthDiff() : 0);
                        cerr << "Mut:          " << mut << endl;
                        cerr << "Off:          " << mut.LengthDiff() << endl;
                        cerr << "Ref:          " << tpl << endl;
                        cerr << "Start,End:    " << start << "," << end << endl;
                        cerr << "PinStart,End: " << pinStart << "," << pinEnd << endl;
                        cerr << "MutStart,End: " << mut.Start() << "," << mut.End() << endl;
                        cerr << "refMutStr:    " << refMutStr << endl;
                        cerr << "c,d:          " << c << "," << d << endl;
                        cerr << "refStr:       " << rStr << endl;
                        cerr << "mutStr:       " << mStr << endl;
                        cerr << "refTpl:       " << ToString(master).substr(c, d - c) << endl;
                        cerr << "mutTpl:       " << ToString(*mTpl) << endl;
                    }

                    // Finally, we should be able to construct a template
                    //  from the mutated string equivalent to Template::Mutate()
                    const int off = rMut->LengthDiff();
                    EXPECT_EQ(end - start + off, mTpl->Length());
                    {
                        Template child(mStr, ModelFactory::Create(mdl, snr));
                        Template applied(rStr, ModelFactory::Create(mdl, snr));
                        applied.ApplyMutations(&rMuts);
                        EXPECT_EQ(child, *mTpl);
                        EXPECT_EQ(child, applied);
                    }
                }
            }
        }
    }
}

TEST(TemplateTest, TestMutatedTemplateEquivalence)
{
#if EXTENSIVE_TESTING
    const int numSamples = 1000;
#else
    const int numSamples = 10;
#endif
    TemplateEquivalence(numSamples, 20, 10);
    TemplateEquivalence(numSamples / 2, 20, 30);
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
        master.ApplyMutation(Mutation::Insertion(len, 'A'));
        EXPECT_EQ(len + 1, master.Length());
        EXPECT_EQ(tpl + A, ToString(master));
        master.ApplyMutation(Mutation::Insertion(0, 'A'));
        EXPECT_EQ(len + 2, master.Length());
        EXPECT_EQ(A + tpl + A, ToString(master));
    }
    // no pinStart but pinEnd
    {
        Template master(tpl, ModelFactory::Create(mdl, snr), 0, len, false, true);
        EXPECT_EQ(len, master.Length());
        EXPECT_EQ(tpl, ToString(master));
        master.ApplyMutation(Mutation::Insertion(0, 'A'));
        EXPECT_EQ(len, master.Length());
        EXPECT_EQ(tpl, ToString(master));
        // the coords are now 1..6, so a new terminal mutation is at len+1
        master.ApplyMutation(Mutation::Insertion(len + 1, 'A'));
        EXPECT_EQ(len + 1, master.Length());
        EXPECT_EQ(tpl + A, ToString(master));
    }
    // pinStart but no pinEnd
    {
        Template master(tpl, ModelFactory::Create(mdl, snr), 0, len, true, false);
        EXPECT_EQ(len, master.Length());
        EXPECT_EQ(tpl, ToString(master));
        master.ApplyMutation(Mutation::Insertion(len, 'A'));
        EXPECT_EQ(len, master.Length());
        EXPECT_EQ(tpl, ToString(master));
        master.ApplyMutation(Mutation::Insertion(0, 'A'));
        EXPECT_EQ(len + 1, master.Length());
        EXPECT_EQ(A + tpl, ToString(master));
    }
    // no pinStart or pinEnd
    {
        Template master(tpl, ModelFactory::Create(mdl, snr), 0, len, false, false);
        EXPECT_EQ(len, master.Length());
        EXPECT_EQ(tpl, ToString(master));
        master.ApplyMutation(Mutation::Insertion(len, 'A'));
        EXPECT_EQ(len, master.Length());
        EXPECT_EQ(tpl, ToString(master));
        master.ApplyMutation(Mutation::Insertion(0, 'A'));
        EXPECT_EQ(len, master.Length());
        EXPECT_EQ(tpl, ToString(master));
    }
}

TEST(TemplateTest, NullTemplate)
{
    const string tpl = "ACGT";
    const size_t len = tpl.length();
    const Mutation mut = Mutation::Deletion(0, 1);

    ASSERT_THROW(Template("A", ModelFactory::Create(mdl, snr), 0, 1, true, true), TemplateTooSmall);

    ASSERT_NO_THROW(Template("AA", ModelFactory::Create(mdl, snr), 0, 2, true, true));
}

TEST(TemplateTest, P6SiteNormalParameters)
{
    const string tpl = "ACGATACATACGATCGA";
    const SNR snr(10, 7, 5, 11);
    auto mdl = "P6-C4";
    Template tester(tpl, ModelFactory::Create(mdl, snr));
    auto results = tester.NormalParameters();

    EXPECT_EQ(-9.3915588824261888, results.first);
    EXPECT_EQ(30.392545575324248, results.second);
}

}  // namespace TemplateTests
