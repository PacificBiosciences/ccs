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

// Author: Lance Hepler

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <string>
#include <tuple>

#include <pacbio/consensus/Integrator.h>
#include <pacbio/consensus/Polish.h>
#include <pacbio/consensus/Read.h>
#include <pacbio/consensus/Sequence.h>

using namespace PacBio::Consensus;  // NOLINT

namespace {

Read MkRead(const std::string& seq, const SNR& snr, const std::string& mdl)
{
    std::vector<uint8_t> cov(0, seq.length());
    return Read("NA", seq, cov, cov, snr, mdl);
}

const SNR snr(10, 7, 5, 11);

TEST(PolishTest, MonoBasic)
{
    MonoMolecularIntegrator ai("GCGTCGT", IntegratorConfig(), snr, "P6-C4");

    ai.AddRead(MappedRead(MkRead("ACGTACGT", snr, "P6-C4"), StrandEnum::FORWARD, 0, 7, true, true));
    ai.AddRead(MappedRead(MkRead("ACGACGT", snr, "P6-C4"), StrandEnum::FORWARD, 0, 7, true, true));
    ai.AddRead(MappedRead(MkRead("ACGACGT", snr, "P6-C4"), StrandEnum::FORWARD, 0, 7, true, true));

    Polish(&ai, PolishConfig());

    EXPECT_EQ("ACGACGT", std::string(ai));
}

TEST(PolishTest, MultiBasic)
{
    MultiMolecularIntegrator ai("GCGTCGT", IntegratorConfig());

    ai.AddRead(MappedRead(MkRead("ACGTACGT", snr, "P6-C4"), StrandEnum::FORWARD, 0, 7, true, true));
    ai.AddRead(MappedRead(MkRead(ReverseComplement("ACGACGT"), snr, "P6-C4"), StrandEnum::REVERSE,
                          0, 7, true, true));
    ai.AddRead(MappedRead(MkRead("ACGACGT", snr, "P6-C4"), StrandEnum::FORWARD, 0, 7, true, true));

    bool polished;
    size_t nTested, nApplied;

    std::tie(polished, nTested, nApplied) = Polish(&ai, PolishConfig());

    EXPECT_TRUE(polished);
    EXPECT_EQ("ACGACGT", std::string(ai));
}
}
