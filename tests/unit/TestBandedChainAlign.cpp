// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
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

// Author: Derek Barnett

#include <gtest/gtest.h>

#include <pacbio/align/BandedChainAlignment.h>
#include <pacbio/align/internal/BCAlignBlocks.h>
#include <pacbio/align/internal/BCAlignImpl.h>

TEST(StandardGlobalAlignBlockTest, Basic)
{
    using Config = PacBio::Align::BandedChainAlignConfig;
    using Alignment = PacBio::Align::BandedChainAlignment;
    using Block = PacBio::Align::Internal::StandardGlobalAlignBlock;

    const Config config = Config::Default();
    Block block{config};

    {  // complete sequence match
        const char* t = "ATT";
        const char* q = "ATT";
        const size_t tLen = 3;
        const size_t qLen = 3;
        const auto cigar = block.Align(t, tLen, q, qLen);
        const auto align = Alignment{config, t, tLen, q, qLen, cigar};
        EXPECT_EQ("3=", cigar.ToStdString());
        EXPECT_EQ(6, align.Score());
    }
    {  // ensure gap at end (query)
        const char* t = "ATT";
        const char* q = "AT";
        const size_t tLen = 3;
        const size_t qLen = 2;
        const auto cigar = block.Align(t, tLen, q, qLen);
        const auto align = Alignment{config, t, tLen, q, qLen, cigar};
        EXPECT_EQ("2=1D", cigar.ToStdString());
        EXPECT_EQ(4, align.Score());  // end-gaps free
    }
    {  // ensure gap at end (target)
        const char* t = "AT";
        const char* q = "ATT";
        const size_t tLen = 2;
        const size_t qLen = 3;
        const auto cigar = block.Align(t, tLen, q, qLen);
        const auto align = Alignment{config, t, tLen, q, qLen, cigar};
        EXPECT_EQ("2=1I", cigar.ToStdString());
        EXPECT_EQ(4, align.Score());  // end-gaps free
    }
    {  // query end insertion
        const char* t = "GA";
        const char* q = "GAT";
        const size_t tLen = 2;
        const size_t qLen = 3;
        const auto cigar = block.Align(t, tLen, q, qLen);
        const auto align = Alignment{config, t, tLen, q, qLen, cigar};
        EXPECT_EQ("2=1I", cigar.ToStdString());
        EXPECT_EQ(4, align.Score());  // end-gaps free
    }
    {  // query end deletion
        const char* t = "GAT";
        const char* q = "GA";
        const size_t tLen = 3;
        const size_t qLen = 2;
        const auto cigar = block.Align(t, tLen, q, qLen);
        const auto align = Alignment{config, t, tLen, q, qLen, cigar};
        EXPECT_EQ("2=1D", cigar.ToStdString());
        EXPECT_EQ(4, align.Score());  // end-gaps free
    }
    {  // query begin insertion
        const char* t = "GA";
        const char* q = "TGA";
        const size_t tLen = 2;
        const size_t qLen = 3;
        const auto cigar = block.Align(t, tLen, q, qLen);
        const auto align = Alignment{config, t, tLen, q, qLen, cigar};
        EXPECT_EQ("1I2=", cigar.ToStdString());
        EXPECT_EQ(4, align.Score());  // end-gaps free
    }
    {  // query begin deletion
        const char* t = "TGA";
        const char* q = "GA";
        const size_t tLen = 3;
        const size_t qLen = 2;
        const auto cigar = block.Align(t, tLen, q, qLen);
        const auto align = Alignment{config, t, tLen, q, qLen, cigar};
        EXPECT_EQ("1D2=", cigar.ToStdString());
        EXPECT_EQ(4, align.Score());  // end-gaps free
    }
    {  // (slightly longer) internal insertion
        const char* t = "GATTACA";
        const char* q = "GATTTACA";
        const size_t tLen = 7;
        const size_t qLen = 8;
        const auto cigar = block.Align(t, tLen, q, qLen);
        const auto align = Alignment{config, t, tLen, q, qLen, cigar};
        EXPECT_EQ("2=1I5=", cigar.ToStdString());
        EXPECT_EQ(12, align.Score());
    }
}

TEST(StandardGlobalAlignBlockTest, LargeGapTest)
{
    const std::string target =
        "AACGATTTTATGATGGCATGTGACATGTATTTCCGTTGGGGGCATTTTAATAAGTGAGGA"
        "AGTGATAGGAAGTGACCAGATAATACATATATGTTCTGTACTCTCTTGCGCATTTTGATT"
        "GTTGACTGAGTAACCAGACAGTTGATGTGCACGATTTCCCCTCGCCCTAACAGACGTGGG"
        "CGGGGGCACCCGGTGTGTGCAATATACGAAAAAAAAGCCCGTACTTTCGTACGAGCTCTT"
        "CTTTAAATATGGCGGTGAGGGGGGGATTGACTCGCTTCGCTCGCCCTGCGGGCAGCCCGC"
        "TCACTGCGTTCACGGTCTGTCCAACTGGCTGTCGCCAGTTGTCGACCCCCGGTCGGGGCT"
        "TCTCATCCCCCCGGTGTGTGCAATACACGAAAAAAAAGCCCGTACTTTCGTACGAGCTCT"
        "TCTTTAAATATGGCGGTGAGGGGGGGATTGACTCGCTTCGCTCGCCCTGCGGGCAGCCCG"
        "CTCACTGCGTTCACGGTCTGTCCAACTGGCTGTCGCCAGTTGTCGAACCCCGGTCGGGGC"
        "TTCTCATCCCCCCGGTGTGTGCAATATACGAAAAAAAAGCCCGTACTTTCGTACGAGCTC"
        "TTCTTTAAATATGGCGGTGAGGGGGGGATTGACTCGCTTCGCTCGCCCTGCGGGCAGCCC"
        "GCTCACTGCGTTCACGGTCTGTCCAACTGGCTGTCGCCAGTTGTCGAACCCCGGTCGGGG"
        "CTTCTCATCCCCCCGGTGTGTGCAATATACGAAAAAAAAGCCCGTACTTTCGTACGAGCT"
        "CTTCTTTAAATATGGCGGTGAGGGGGGGATTCGAACCCCCGATACGTTGCCGTATACACA"
        "CTTTCCAGGCGTGCTCCTTCAGCCACTCGGACACCTCACCAAATTGTCGTTCCTGTCTTG"
        "CTGGAACGGGCGCTAATTTAGGGAAATCATGACCTGAGGTCAACAAACTTTTTGAAAAAA"
        "TCGCGCGTTTATTCAAACTTCAATCAATGTGTGGTTTTAATAAGCGAAAT";

    const std::string query =
        "AACGATTTTATGATGGCATGTGACATGTATTTCCGTTGGGGGCATTTTAATAAGTGAGGA"
        "AGTGATAGGAAGTGACCAGATAATACATATATGTTCTGTACTCTCTTGCGCATTTTGATT"
        "GTTGACTGAGTAACCAGACAGTTGATGTGCACGATTTCCCCTCGCCCTAACAGACGTGGG"
        "CGGGGGCACCCGGTGTGTGCAATATACGAAAAAAAAGCCCGTACTTTCGTACGAGCTCTT"
        "CTTTAAATATGGCGGTGAGGGGGGGATTGACTCGCTTCGCTCGCCCTGCGGGCAGCCCGC"
        "TCACTGCGTTCACGGTCTGTCCAACTGGCTGTCGCCAGTTGTCGACCCCCGGTCGGGGCT"
        "TCTCATCCCCCCGGTGTGTGCAATACACGAAAAAAAAGCCCGTACTTTCGTACGAGCTCT"
        "TCTTTAAATATGGCGGTGAGGGGGGGATTGACTCGCTTCGCTCGCCCTGCGGGCAGCCCG"
        "CTCACTGCGTTCACGGTCTGTCCAACTGGCTGTCGCCAGTTGTCGAACCCCGGTCGGGGC"
        "TTCTCATCCCCCCGGTGTGTGCAATATACGAAAAAAAAGCCCGTACTTTCGTACGAGCTC"
        "TTCTTTAAATATGGCGGTGAGGGGGGGATTGACTCGCTTCGCTCGCCCTGCGGGCAGCCC"
        "GCTCACTGCGTTCACGGTCTGTCCAACTGGCTGTCGCCAGTTGTCGAACCCCGGTCGGGG"
        "CTTCTCATCCCCCCGGTGTGTGCAATATACGAAAAAAAAGCCCGTACTTTCGTACGAGCT"
        "CTTCTTTAAATATGGCGGTGAGGGGGGGATTGACTCGCTTCGCTCGCCCTGCGGGCAGCC"
        "CGCTCACTGCGTTCACGGTCTGTCCAACTGGCTGTCGCCAGTTGTCGAACCCCGGTCGGG"
        "GCTTCTCATCCCCCCGGTGTGTGCAATATACGAAAAAAAAGCCCGTACTTTCGTACGAGC"
        "TCTTCTTTAAATATGGCGGTGAGGGGGGGATTCGAACCCCCGATACGTTGCCGTATACAC"
        "ACTTTCCAGGCGTGCTCCTTCAGCCACTCGGACACCTCACCAAATTGTCGTTCCTGTCTT"
        "GCTGGAACGGGCGCTAATTTAGGGAAATCATGACCTGAGGTCAACAAACTTTTTGAAAAA"
        "ATCGCGCGTTTATTCAAACTTCAATCAATGTGTGGTTTTAATAAGCGAAAT";

    using Config = PacBio::Align::BandedChainAlignConfig;
    using Block = PacBio::Align::Internal::StandardGlobalAlignBlock;
    const Config config{2, -1, -2, -1, 1};
    Block block{config};
    const auto cigar = block.Align(target.c_str(), target.size(), query.c_str(), query.size());
    EXPECT_EQ("386=181I624=", cigar.ToStdString());
}

TEST(BandedGlobalAlignBlockTest, Basic)
{
    using Config = PacBio::Align::BandedChainAlignConfig;
    using Alignment = PacBio::Align::BandedChainAlignment;
    using Block = PacBio::Align::Internal::BandedGlobalAlignBlock;

    const char* t = "GATTACAT";
    const char* q = "GATTACAT";
    const size_t tLen = 8;
    const size_t qLen = 8;
    const auto seed = PacBio::Align::Seed{0, 0, 8};
    Config config = Config::Default();
    config.bandExtend_ = 2;
    Block block{config};
    const auto cigar = block.Align(t, q, seed);
    const auto align = Alignment{config, t, tLen, q, qLen, cigar};

    EXPECT_EQ("8=", cigar.ToStdString());
    EXPECT_EQ(16, align.Score());
}

TEST(BandedGlobalAlignBlockTest, Align)
{
    using Config = PacBio::Align::BandedChainAlignConfig;
    using Alignment = PacBio::Align::BandedChainAlignment;
    using Block = PacBio::Align::Internal::BandedGlobalAlignBlock;
    using Seed = PacBio::Align::Seed;

    Config config = Config::Default();
    config.bandExtend_ = 2;
    Block block{config};

    {
        const char* t = "ATAGAT";
        const char* q = "ATGT";
        const size_t tLen = 6;
        const size_t qLen = 4;
        const Seed seed{0, 0, 6, 4};

        const auto cigar = block.Align(t, q, seed);
        const auto align = Alignment{config, t, tLen, q, qLen, cigar};

        // ATAGAT
        // AT-G-T

        EXPECT_EQ("2=1D1=1D1=", cigar.ToStdString());
        EXPECT_EQ(4, align.Score());
    }
    {
        const char* t = "AAAAAATTTTTGGG";
        const char* q = "TTTTTTTTGGGGGGGG";
        const size_t tLen = 14;
        const size_t qLen = 16;
        const Seed seed{0, 0, 14, 16};  // no offset

        const auto cigar = block.Align(t, q, seed);
        const auto align = Alignment{config, t, tLen, q, qLen, cigar};

        // Expected:
        // 2D4X4=1X3=4I
        // AAAAAATTTTTGGG----
        // --TTTTTTTTGGGGGGGG

        EXPECT_EQ("2D4X4=1X3=4I", cigar.ToStdString());
        EXPECT_EQ(9, align.Score());  // end-gaps free
    }
}

TEST(BandedAligner, CigarStitching)
{
    using Cigar = PacBio::Data::Cigar;
    using Config = PacBio::Align::BandedChainAlignConfig;
    using Aligner = PacBio::Align::Internal::BandedChainAlignerImpl;

    const Config config = Config::Default();

    {  // simple
        Cigar global("5=");
        Cigar local("2=");

        Aligner a{config};
        a.StitchCigars(&global, std::move(local));

        EXPECT_EQ("7=", global.ToStdString());
    }
    {  // different at edge
        Cigar global("5=3D");
        Cigar local("2=1X");

        Aligner a{config};
        a.StitchCigars(&global, std::move(local));

        EXPECT_EQ("5=3D2=1X", global.ToStdString());
    }
}

TEST(BandedAligner, AlignSeeds)
{
    using Config = PacBio::Align::BandedChainAlignConfig;
    using Seed = PacBio::Align::Seed;

    {
        Config config = Config::Default();
        config.bandExtend_ = 2;
        const std::string target = "CGAATCCATCCCACACA";
        const std::string query = "GGCGATNNNCATGGCACA";
        const auto seeds =
            std::vector<Seed>{Seed{0, 2, 5, 6}, Seed{6, 9, 9, 12}, Seed{11, 14, 17, 16}};

        const auto result = PacBio::Align::BandedChainAlign(target, query, seeds, config);

        EXPECT_EQ("--CGAATC--CATCCCACACA", result.alignedTarget_);
        EXPECT_EQ("GGCG-ATNNNCATGGCACA--", result.alignedQuery_);
        EXPECT_EQ("2I2=1D2=1X2I3=2X4=2D", result.cigar_.ToStdString());
        EXPECT_EQ(14, result.Score());  // end-gaps free
    }
}
