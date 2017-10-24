// Copyright (c) 2011-2017, Pacific Biosciences of California, Inc.
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

// Author: Lance Hepler

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <iostream>
#include <string>
#include <tuple>

using std::string;

#include <pacbio/consensus/Integrator.h>
#include <pacbio/consensus/Polish.h>
#include <pacbio/data/Read.h>
#include <pacbio/data/Sequence.h>

using namespace PacBio::Consensus;  // NOLINT
using namespace PacBio::Data;       // NOLINT

namespace PolishTests {

Read MkRead(const std::string& seq, const SNR& snr, const std::string& mdl)
{
    std::vector<uint8_t> cov(seq.length(), 0);
    return Read("NA", seq, cov, cov, snr, mdl);
}

const SNR snr(10, 7, 5, 11);
const string mdl("P6-C4");

TEST(PolishTest, Basic)
{
    Integrator ai("GCGTCGT", IntegratorConfig());

    ai.AddRead(MappedRead(MkRead("ACGTACGT", snr, mdl), StrandType::FORWARD, 0, 7, true, true));
    ai.AddRead(MappedRead(MkRead(ReverseComplement("ACGACGT"), snr, mdl), StrandType::REVERSE, 0, 7,
                          true, true));
    ai.AddRead(MappedRead(MkRead("ACGACGT", snr, mdl), StrandType::FORWARD, 0, 7, true, true));

    const auto result = Polish(&ai, PolishConfig());

    EXPECT_TRUE(result.hasConverged);
    EXPECT_EQ("ACGACGT", string(ai));
}

TEST(PolishTest, DiTriRepeat)
{
    //                       1  2  31 2 3
    const string tpl = "ACGTCAGCAGCAGAGAGTGCA";
    //                       1  2  3  41 2 3 4
    const string read = "ACGTCAGCAGCAGCAGAGAGAGTGCA";
    Integrator ai(tpl, IntegratorConfig());

    ai.AddRead(
        MappedRead(MkRead(read, snr, mdl), StrandType::FORWARD, 0, tpl.length(), true, true));
    ai.AddRead(MappedRead(MkRead(ReverseComplement(read), snr, mdl), StrandType::REVERSE, 0,
                          tpl.length(), true, true));
    ai.AddRead(
        MappedRead(MkRead(read, snr, mdl), StrandType::FORWARD, 0, tpl.length(), true, true));
    // ai.AddRead(MappedRead(MkRead(ReverseComplement(tpl), snr, mdl), StrandType::REVERSE, 0, tpl.length(), true, true));
    // ai.AddRead(MappedRead(MkRead(tpl, snr, mdl), StrandType::FORWARD, 0, tpl.length(), true, true));

    if (false) {
        const auto noresult = Polish(&ai, PolishConfig());

        EXPECT_TRUE(noresult.hasConverged);
        EXPECT_EQ(tpl, string(ai));
    }

    const auto result = PolishRepeats(&ai, RepeatConfig());

    EXPECT_TRUE(result.hasConverged);
    EXPECT_EQ(read, string(ai));
}
}  // namespace PolishTests
