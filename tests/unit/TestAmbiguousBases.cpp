// Author: David Seifert

#include <gtest/gtest.h>

#include <pacbio/data/internal/BaseEncoding.h>

using namespace PacBio::Data::detail;

TEST(AmbiguousBasesTest, TestNCBI2naValid)
{
    const auto testBaseA = NCBI2na::FromASCII('A');
    const auto testBaseC = NCBI2na::FromASCII('C');
    const auto testBaseG = NCBI2na::FromASCII('G');
    const auto testBaseT = NCBI2na::FromASCII('T');

    EXPECT_EQ(testBaseA.IsValid(), true);
    EXPECT_EQ(testBaseC.IsValid(), true);
    EXPECT_EQ(testBaseG.IsValid(), true);
    EXPECT_EQ(testBaseT.IsValid(), true);
}

TEST(AmbiguousBasesTest, TestNCBI4naValid)
{
    const auto testBaseA = NCBI4na::FromASCII('A');
    const auto testBaseC = NCBI4na::FromASCII('C');
    const auto testBaseM = NCBI4na::FromASCII('M');
    const auto testBaseG = NCBI4na::FromASCII('G');
    const auto testBaseR = NCBI4na::FromASCII('R');
    const auto testBaseS = NCBI4na::FromASCII('S');
    const auto testBaseV = NCBI4na::FromASCII('V');
    const auto testBaseT = NCBI4na::FromASCII('T');
    const auto testBaseW = NCBI4na::FromASCII('W');
    const auto testBaseY = NCBI4na::FromASCII('Y');
    const auto testBaseH = NCBI4na::FromASCII('H');
    const auto testBaseK = NCBI4na::FromASCII('K');
    const auto testBaseD = NCBI4na::FromASCII('D');
    const auto testBaseB = NCBI4na::FromASCII('B');
    const auto testBaseN = NCBI4na::FromASCII('N');

    EXPECT_EQ(testBaseA.IsValid(), true);
    EXPECT_EQ(testBaseC.IsValid(), true);
    EXPECT_EQ(testBaseM.IsValid(), true);
    EXPECT_EQ(testBaseG.IsValid(), true);
    EXPECT_EQ(testBaseR.IsValid(), true);
    EXPECT_EQ(testBaseS.IsValid(), true);
    EXPECT_EQ(testBaseV.IsValid(), true);
    EXPECT_EQ(testBaseT.IsValid(), true);
    EXPECT_EQ(testBaseW.IsValid(), true);
    EXPECT_EQ(testBaseY.IsValid(), true);
    EXPECT_EQ(testBaseH.IsValid(), true);
    EXPECT_EQ(testBaseK.IsValid(), true);
    EXPECT_EQ(testBaseD.IsValid(), true);
    EXPECT_EQ(testBaseB.IsValid(), true);
    EXPECT_EQ(testBaseN.IsValid(), true);
}

TEST(AmbiguousBasesTest, TestNCBI2naInvalidAmbig)
{
    // check ambiguous bases
    const auto testBaseM = NCBI2na::FromASCII('M');
    const auto testBaseR = NCBI2na::FromASCII('R');
    const auto testBaseS = NCBI2na::FromASCII('S');
    const auto testBaseV = NCBI2na::FromASCII('V');
    const auto testBaseW = NCBI2na::FromASCII('W');
    const auto testBaseY = NCBI2na::FromASCII('Y');
    const auto testBaseH = NCBI2na::FromASCII('H');
    const auto testBaseK = NCBI2na::FromASCII('K');
    const auto testBaseD = NCBI2na::FromASCII('D');
    const auto testBaseB = NCBI2na::FromASCII('B');
    const auto testBaseN = NCBI2na::FromASCII('N');

    EXPECT_EQ(testBaseM.IsValid(), false);
    EXPECT_EQ(testBaseR.IsValid(), false);
    EXPECT_EQ(testBaseS.IsValid(), false);
    EXPECT_EQ(testBaseV.IsValid(), false);
    EXPECT_EQ(testBaseW.IsValid(), false);
    EXPECT_EQ(testBaseY.IsValid(), false);
    EXPECT_EQ(testBaseH.IsValid(), false);
    EXPECT_EQ(testBaseK.IsValid(), false);
    EXPECT_EQ(testBaseD.IsValid(), false);
    EXPECT_EQ(testBaseB.IsValid(), false);
    EXPECT_EQ(testBaseN.IsValid(), false);
}

TEST(AmbiguousBasesTest, TestNCBI2naInvalidGarbage)
{
    // check absolutely garbage bases
    const auto testBase0 = NCBI2na::FromASCII('\0');
    const auto testBaseQ = NCBI2na::FromASCII('Q');
    const auto testBaseX = NCBI2na::FromASCII('X');
    const auto testBaseZ = NCBI2na::FromASCII('Z');

    EXPECT_EQ(testBase0.IsValid(), false);
    EXPECT_EQ(testBaseQ.IsValid(), false);
    EXPECT_EQ(testBaseX.IsValid(), false);
    EXPECT_EQ(testBaseZ.IsValid(), false);
}

TEST(AmbiguousBasesTest, TestNCBI4naInvalidGarbage)
{
    // check absolutely garbage bases
    EXPECT_ANY_THROW({ const auto testBase0 = NCBI4na::FromASCII('\0'); });
    EXPECT_ANY_THROW({ const auto testBaseQ = NCBI4na::FromASCII('Q'); });
    EXPECT_ANY_THROW({ const auto testBaseX = NCBI4na::FromASCII('X'); });
    EXPECT_ANY_THROW({ const auto testBaseZ = NCBI4na::FromASCII('Z'); });
}

TEST(AmbiguousBasesTest, TestNCBI2naASCII)
{
    const auto testBaseA = NCBI2na::FromASCII('A');
    const auto testBaseC = NCBI2na::FromASCII('C');
    const auto testBaseG = NCBI2na::FromASCII('G');
    const auto testBaseT = NCBI2na::FromASCII('T');

    EXPECT_EQ(testBaseA.GetASCII(), 'A');
    EXPECT_EQ(testBaseC.GetASCII(), 'C');
    EXPECT_EQ(testBaseG.GetASCII(), 'G');
    EXPECT_EQ(testBaseT.GetASCII(), 'T');
}

TEST(AmbiguousBasesTest, TestNCBI4naASCII)
{
    const auto testBaseA = NCBI4na::FromASCII('A');
    const auto testBaseC = NCBI4na::FromASCII('C');
    const auto testBaseM = NCBI4na::FromASCII('M');
    const auto testBaseG = NCBI4na::FromASCII('G');
    const auto testBaseR = NCBI4na::FromASCII('R');
    const auto testBaseS = NCBI4na::FromASCII('S');
    const auto testBaseV = NCBI4na::FromASCII('V');
    const auto testBaseT = NCBI4na::FromASCII('T');
    const auto testBaseW = NCBI4na::FromASCII('W');
    const auto testBaseY = NCBI4na::FromASCII('Y');
    const auto testBaseH = NCBI4na::FromASCII('H');
    const auto testBaseK = NCBI4na::FromASCII('K');
    const auto testBaseD = NCBI4na::FromASCII('D');
    const auto testBaseB = NCBI4na::FromASCII('B');
    const auto testBaseN = NCBI4na::FromASCII('N');

    EXPECT_EQ(testBaseA.GetASCII(), 'A');
    EXPECT_EQ(testBaseC.GetASCII(), 'C');
    EXPECT_EQ(testBaseM.GetASCII(), 'M');
    EXPECT_EQ(testBaseG.GetASCII(), 'G');
    EXPECT_EQ(testBaseR.GetASCII(), 'R');
    EXPECT_EQ(testBaseS.GetASCII(), 'S');
    EXPECT_EQ(testBaseV.GetASCII(), 'V');
    EXPECT_EQ(testBaseT.GetASCII(), 'T');
    EXPECT_EQ(testBaseW.GetASCII(), 'W');
    EXPECT_EQ(testBaseY.GetASCII(), 'Y');
    EXPECT_EQ(testBaseH.GetASCII(), 'H');
    EXPECT_EQ(testBaseK.GetASCII(), 'K');
    EXPECT_EQ(testBaseD.GetASCII(), 'D');
    EXPECT_EQ(testBaseB.GetASCII(), 'B');
    EXPECT_EQ(testBaseN.GetASCII(), 'N');
}

TEST(AmbiguousBasesTest, TestNCBI2naToNCBI4naBijection)
{
    const auto testBaseNCBI2naA = NCBI2na::FromASCII('A');
    const auto testBaseNCBI2naC = NCBI2na::FromASCII('C');
    const auto testBaseNCBI2naG = NCBI2na::FromASCII('G');
    const auto testBaseNCBI2naT = NCBI2na::FromASCII('T');

    const auto testBaseNCBI4naA = testBaseNCBI2naA.GetNCBI4na();
    const auto testBaseNCBI4naC = testBaseNCBI2naC.GetNCBI4na();
    const auto testBaseNCBI4naG = testBaseNCBI2naG.GetNCBI4na();
    const auto testBaseNCBI4naT = testBaseNCBI2naT.GetNCBI4na();

    // 1. Test equal
    EXPECT_EQ(testBaseNCBI2naA.IsEqual(testBaseNCBI4naA.GetNCBI2na()), true);
    EXPECT_EQ(testBaseNCBI2naC.IsEqual(testBaseNCBI4naC.GetNCBI2na()), true);
    EXPECT_EQ(testBaseNCBI2naG.IsEqual(testBaseNCBI4naG.GetNCBI2na()), true);
    EXPECT_EQ(testBaseNCBI2naT.IsEqual(testBaseNCBI4naT.GetNCBI2na()), true);

    // 2. Test unequal
    // A
    EXPECT_EQ(testBaseNCBI2naA.IsEqual(testBaseNCBI4naC.GetNCBI2na()), false);
    EXPECT_EQ(testBaseNCBI2naA.IsEqual(testBaseNCBI4naG.GetNCBI2na()), false);
    EXPECT_EQ(testBaseNCBI2naA.IsEqual(testBaseNCBI4naT.GetNCBI2na()), false);

    // C
    EXPECT_EQ(testBaseNCBI2naC.IsEqual(testBaseNCBI4naA.GetNCBI2na()), false);
    EXPECT_EQ(testBaseNCBI2naC.IsEqual(testBaseNCBI4naG.GetNCBI2na()), false);
    EXPECT_EQ(testBaseNCBI2naC.IsEqual(testBaseNCBI4naT.GetNCBI2na()), false);

    // G
    EXPECT_EQ(testBaseNCBI2naG.IsEqual(testBaseNCBI4naA.GetNCBI2na()), false);
    EXPECT_EQ(testBaseNCBI2naG.IsEqual(testBaseNCBI4naC.GetNCBI2na()), false);
    EXPECT_EQ(testBaseNCBI2naG.IsEqual(testBaseNCBI4naT.GetNCBI2na()), false);

    // T
    EXPECT_EQ(testBaseNCBI2naT.IsEqual(testBaseNCBI4naA.GetNCBI2na()), false);
    EXPECT_EQ(testBaseNCBI2naT.IsEqual(testBaseNCBI4naC.GetNCBI2na()), false);
    EXPECT_EQ(testBaseNCBI2naT.IsEqual(testBaseNCBI4naG.GetNCBI2na()), false);
}

TEST(AmbiguousBasesTest, TestNCBI4naIsPure)
{
    const auto testBaseA = NCBI4na::FromASCII('A');
    const auto testBaseC = NCBI4na::FromASCII('C');
    const auto testBaseM = NCBI4na::FromASCII('M');
    const auto testBaseG = NCBI4na::FromASCII('G');
    const auto testBaseR = NCBI4na::FromASCII('R');
    const auto testBaseS = NCBI4na::FromASCII('S');
    const auto testBaseV = NCBI4na::FromASCII('V');
    const auto testBaseT = NCBI4na::FromASCII('T');
    const auto testBaseW = NCBI4na::FromASCII('W');
    const auto testBaseY = NCBI4na::FromASCII('Y');
    const auto testBaseH = NCBI4na::FromASCII('H');
    const auto testBaseK = NCBI4na::FromASCII('K');
    const auto testBaseD = NCBI4na::FromASCII('D');
    const auto testBaseB = NCBI4na::FromASCII('B');
    const auto testBaseN = NCBI4na::FromASCII('N');

    EXPECT_EQ(testBaseA.IsPure(), true);
    EXPECT_EQ(testBaseC.IsPure(), true);
    EXPECT_EQ(testBaseM.IsPure(), false);
    EXPECT_EQ(testBaseG.IsPure(), true);
    EXPECT_EQ(testBaseR.IsPure(), false);
    EXPECT_EQ(testBaseS.IsPure(), false);
    EXPECT_EQ(testBaseV.IsPure(), false);
    EXPECT_EQ(testBaseT.IsPure(), true);
    EXPECT_EQ(testBaseW.IsPure(), false);
    EXPECT_EQ(testBaseY.IsPure(), false);
    EXPECT_EQ(testBaseH.IsPure(), false);
    EXPECT_EQ(testBaseK.IsPure(), false);
    EXPECT_EQ(testBaseD.IsPure(), false);
    EXPECT_EQ(testBaseB.IsPure(), false);
    EXPECT_EQ(testBaseN.IsPure(), false);
}

TEST(AmbiguousBasesTest, TestNCBI4naIsAmbig)
{
    const auto testBaseA = NCBI4na::FromASCII('A');
    const auto testBaseC = NCBI4na::FromASCII('C');
    const auto testBaseM = NCBI4na::FromASCII('M');
    const auto testBaseG = NCBI4na::FromASCII('G');
    const auto testBaseR = NCBI4na::FromASCII('R');
    const auto testBaseS = NCBI4na::FromASCII('S');
    const auto testBaseV = NCBI4na::FromASCII('V');
    const auto testBaseT = NCBI4na::FromASCII('T');
    const auto testBaseW = NCBI4na::FromASCII('W');
    const auto testBaseY = NCBI4na::FromASCII('Y');
    const auto testBaseH = NCBI4na::FromASCII('H');
    const auto testBaseK = NCBI4na::FromASCII('K');
    const auto testBaseD = NCBI4na::FromASCII('D');
    const auto testBaseB = NCBI4na::FromASCII('B');
    const auto testBaseN = NCBI4na::FromASCII('N');

    EXPECT_EQ(testBaseA.IsAmbig(), false);
    EXPECT_EQ(testBaseC.IsAmbig(), false);
    EXPECT_EQ(testBaseM.IsAmbig(), true);
    EXPECT_EQ(testBaseG.IsAmbig(), false);
    EXPECT_EQ(testBaseR.IsAmbig(), true);
    EXPECT_EQ(testBaseS.IsAmbig(), true);
    EXPECT_EQ(testBaseV.IsAmbig(), true);
    EXPECT_EQ(testBaseT.IsAmbig(), false);
    EXPECT_EQ(testBaseW.IsAmbig(), true);
    EXPECT_EQ(testBaseY.IsAmbig(), true);
    EXPECT_EQ(testBaseH.IsAmbig(), true);
    EXPECT_EQ(testBaseK.IsAmbig(), true);
    EXPECT_EQ(testBaseD.IsAmbig(), true);
    EXPECT_EQ(testBaseB.IsAmbig(), true);
    EXPECT_EQ(testBaseN.IsAmbig(), true);
}

TEST(AmbiguousBasesTest, TestNCBI4naNumContainedBases)
{
    const auto testBaseA = NCBI4na::FromASCII('A');
    const auto testBaseC = NCBI4na::FromASCII('C');
    const auto testBaseM = NCBI4na::FromASCII('M');
    const auto testBaseG = NCBI4na::FromASCII('G');
    const auto testBaseR = NCBI4na::FromASCII('R');
    const auto testBaseS = NCBI4na::FromASCII('S');
    const auto testBaseV = NCBI4na::FromASCII('V');
    const auto testBaseT = NCBI4na::FromASCII('T');
    const auto testBaseW = NCBI4na::FromASCII('W');
    const auto testBaseY = NCBI4na::FromASCII('Y');
    const auto testBaseH = NCBI4na::FromASCII('H');
    const auto testBaseK = NCBI4na::FromASCII('K');
    const auto testBaseD = NCBI4na::FromASCII('D');
    const auto testBaseB = NCBI4na::FromASCII('B');
    const auto testBaseN = NCBI4na::FromASCII('N');

    // haploid
    EXPECT_EQ(testBaseA.NumContainedBases(), 1);
    EXPECT_EQ(testBaseC.NumContainedBases(), 1);
    EXPECT_EQ(testBaseG.NumContainedBases(), 1);
    EXPECT_EQ(testBaseT.NumContainedBases(), 1);

    // diploid
    EXPECT_EQ(testBaseM.NumContainedBases(), 2);
    EXPECT_EQ(testBaseR.NumContainedBases(), 2);
    EXPECT_EQ(testBaseS.NumContainedBases(), 2);
    EXPECT_EQ(testBaseW.NumContainedBases(), 2);
    EXPECT_EQ(testBaseY.NumContainedBases(), 2);
    EXPECT_EQ(testBaseK.NumContainedBases(), 2);

    // triploid
    EXPECT_EQ(testBaseV.NumContainedBases(), 3);
    EXPECT_EQ(testBaseH.NumContainedBases(), 3);
    EXPECT_EQ(testBaseD.NumContainedBases(), 3);
    EXPECT_EQ(testBaseB.NumContainedBases(), 3);

    // tetraploid
    EXPECT_EQ(testBaseN.NumContainedBases(), 4);
}

// 1. perform all tests where intersection is NOT empty
//    i.e., where the NCBI2na base is contained within
//    the NCBI4na base.
TEST(AmbiguousBasesTest, TestNCBI2naContainedInNCBI4na)
{
    // complete NCBI2na space
    const auto testBaseNCBI2naA = NCBI2na::FromASCII('A');
    const auto testBaseNCBI2naC = NCBI2na::FromASCII('C');
    const auto testBaseNCBI2naG = NCBI2na::FromASCII('G');
    const auto testBaseNCBI2naT = NCBI2na::FromASCII('T');

    // complete NCBI4na space
    const auto testBaseNCBI4naA = NCBI4na::FromASCII('A');
    const auto testBaseNCBI4naC = NCBI4na::FromASCII('C');
    const auto testBaseNCBI4naM = NCBI4na::FromASCII('M');
    const auto testBaseNCBI4naG = NCBI4na::FromASCII('G');
    const auto testBaseNCBI4naR = NCBI4na::FromASCII('R');
    const auto testBaseNCBI4naS = NCBI4na::FromASCII('S');
    const auto testBaseNCBI4naV = NCBI4na::FromASCII('V');
    const auto testBaseNCBI4naT = NCBI4na::FromASCII('T');
    const auto testBaseNCBI4naW = NCBI4na::FromASCII('W');
    const auto testBaseNCBI4naY = NCBI4na::FromASCII('Y');
    const auto testBaseNCBI4naH = NCBI4na::FromASCII('H');
    const auto testBaseNCBI4naK = NCBI4na::FromASCII('K');
    const auto testBaseNCBI4naD = NCBI4na::FromASCII('D');
    const auto testBaseNCBI4naB = NCBI4na::FromASCII('B');
    const auto testBaseNCBI4naN = NCBI4na::FromASCII('N');

    // A \in {A, R, W, M, D, V, H, N}
    EXPECT_EQ(testBaseNCBI4naA.Contains(testBaseNCBI2naA), true);
    EXPECT_EQ(testBaseNCBI4naR.Contains(testBaseNCBI2naA), true);
    EXPECT_EQ(testBaseNCBI4naW.Contains(testBaseNCBI2naA), true);
    EXPECT_EQ(testBaseNCBI4naM.Contains(testBaseNCBI2naA), true);
    EXPECT_EQ(testBaseNCBI4naD.Contains(testBaseNCBI2naA), true);
    EXPECT_EQ(testBaseNCBI4naV.Contains(testBaseNCBI2naA), true);
    EXPECT_EQ(testBaseNCBI4naH.Contains(testBaseNCBI2naA), true);
    EXPECT_EQ(testBaseNCBI4naN.Contains(testBaseNCBI2naA), true);

    // C \in {C, Y, S, M, V, H, B, N}
    EXPECT_EQ(testBaseNCBI4naC.Contains(testBaseNCBI2naC), true);
    EXPECT_EQ(testBaseNCBI4naY.Contains(testBaseNCBI2naC), true);
    EXPECT_EQ(testBaseNCBI4naS.Contains(testBaseNCBI2naC), true);
    EXPECT_EQ(testBaseNCBI4naM.Contains(testBaseNCBI2naC), true);
    EXPECT_EQ(testBaseNCBI4naV.Contains(testBaseNCBI2naC), true);
    EXPECT_EQ(testBaseNCBI4naH.Contains(testBaseNCBI2naC), true);
    EXPECT_EQ(testBaseNCBI4naB.Contains(testBaseNCBI2naC), true);
    EXPECT_EQ(testBaseNCBI4naN.Contains(testBaseNCBI2naC), true);

    // G \in {G, R, S, K, D, V, B, N}
    EXPECT_EQ(testBaseNCBI4naG.Contains(testBaseNCBI2naG), true);
    EXPECT_EQ(testBaseNCBI4naR.Contains(testBaseNCBI2naG), true);
    EXPECT_EQ(testBaseNCBI4naS.Contains(testBaseNCBI2naG), true);
    EXPECT_EQ(testBaseNCBI4naK.Contains(testBaseNCBI2naG), true);
    EXPECT_EQ(testBaseNCBI4naD.Contains(testBaseNCBI2naG), true);
    EXPECT_EQ(testBaseNCBI4naV.Contains(testBaseNCBI2naG), true);
    EXPECT_EQ(testBaseNCBI4naB.Contains(testBaseNCBI2naG), true);
    EXPECT_EQ(testBaseNCBI4naN.Contains(testBaseNCBI2naG), true);

    // T \in {T, Y, W, K, D, H, B, N}
    EXPECT_EQ(testBaseNCBI4naT.Contains(testBaseNCBI2naT), true);
    EXPECT_EQ(testBaseNCBI4naY.Contains(testBaseNCBI2naT), true);
    EXPECT_EQ(testBaseNCBI4naW.Contains(testBaseNCBI2naT), true);
    EXPECT_EQ(testBaseNCBI4naK.Contains(testBaseNCBI2naT), true);
    EXPECT_EQ(testBaseNCBI4naD.Contains(testBaseNCBI2naT), true);
    EXPECT_EQ(testBaseNCBI4naH.Contains(testBaseNCBI2naT), true);
    EXPECT_EQ(testBaseNCBI4naB.Contains(testBaseNCBI2naT), true);
    EXPECT_EQ(testBaseNCBI4naN.Contains(testBaseNCBI2naT), true);
}

// 2. perform all tests where intersection is empty
//    i.e., where the NCBI2na base is NOT contained
//    within the NCBI4na base.
TEST(AmbiguousBasesTest, TestNCBI2naNotContainedInNCBI4na)
{
    // complete NCBI2na space
    const auto testBaseNCBI2naA = NCBI2na::FromASCII('A');
    const auto testBaseNCBI2naC = NCBI2na::FromASCII('C');
    const auto testBaseNCBI2naG = NCBI2na::FromASCII('G');
    const auto testBaseNCBI2naT = NCBI2na::FromASCII('T');

    // complete NCBI4na space
    const auto testBaseNCBI4naA = NCBI4na::FromASCII('A');
    const auto testBaseNCBI4naC = NCBI4na::FromASCII('C');
    const auto testBaseNCBI4naM = NCBI4na::FromASCII('M');
    const auto testBaseNCBI4naG = NCBI4na::FromASCII('G');
    const auto testBaseNCBI4naR = NCBI4na::FromASCII('R');
    const auto testBaseNCBI4naS = NCBI4na::FromASCII('S');
    const auto testBaseNCBI4naV = NCBI4na::FromASCII('V');
    const auto testBaseNCBI4naT = NCBI4na::FromASCII('T');
    const auto testBaseNCBI4naW = NCBI4na::FromASCII('W');
    const auto testBaseNCBI4naY = NCBI4na::FromASCII('Y');
    const auto testBaseNCBI4naH = NCBI4na::FromASCII('H');
    const auto testBaseNCBI4naK = NCBI4na::FromASCII('K');
    const auto testBaseNCBI4naD = NCBI4na::FromASCII('D');
    const auto testBaseNCBI4naB = NCBI4na::FromASCII('B');
    const auto testBaseNCBI4naN = NCBI4na::FromASCII('N');

    // A \not\in {C, G, T, Y, S, K, B}
    EXPECT_EQ(testBaseNCBI4naC.Contains(testBaseNCBI2naA), false);
    EXPECT_EQ(testBaseNCBI4naG.Contains(testBaseNCBI2naA), false);
    EXPECT_EQ(testBaseNCBI4naT.Contains(testBaseNCBI2naA), false);
    EXPECT_EQ(testBaseNCBI4naY.Contains(testBaseNCBI2naA), false);
    EXPECT_EQ(testBaseNCBI4naS.Contains(testBaseNCBI2naA), false);
    EXPECT_EQ(testBaseNCBI4naK.Contains(testBaseNCBI2naA), false);
    EXPECT_EQ(testBaseNCBI4naB.Contains(testBaseNCBI2naA), false);

    // C \not\in {A, G, T, R, W, K, D}
    EXPECT_EQ(testBaseNCBI4naA.Contains(testBaseNCBI2naC), false);
    EXPECT_EQ(testBaseNCBI4naG.Contains(testBaseNCBI2naC), false);
    EXPECT_EQ(testBaseNCBI4naT.Contains(testBaseNCBI2naC), false);
    EXPECT_EQ(testBaseNCBI4naR.Contains(testBaseNCBI2naC), false);
    EXPECT_EQ(testBaseNCBI4naW.Contains(testBaseNCBI2naC), false);
    EXPECT_EQ(testBaseNCBI4naK.Contains(testBaseNCBI2naC), false);
    EXPECT_EQ(testBaseNCBI4naD.Contains(testBaseNCBI2naC), false);

    // G \not\in {A, C, T, Y, W, M, H}
    EXPECT_EQ(testBaseNCBI4naA.Contains(testBaseNCBI2naG), false);
    EXPECT_EQ(testBaseNCBI4naC.Contains(testBaseNCBI2naG), false);
    EXPECT_EQ(testBaseNCBI4naT.Contains(testBaseNCBI2naG), false);
    EXPECT_EQ(testBaseNCBI4naY.Contains(testBaseNCBI2naG), false);
    EXPECT_EQ(testBaseNCBI4naW.Contains(testBaseNCBI2naG), false);
    EXPECT_EQ(testBaseNCBI4naM.Contains(testBaseNCBI2naG), false);
    EXPECT_EQ(testBaseNCBI4naH.Contains(testBaseNCBI2naG), false);

    // T \not\in {A, C, G, R, S, M, V}
    EXPECT_EQ(testBaseNCBI4naA.Contains(testBaseNCBI2naT), false);
    EXPECT_EQ(testBaseNCBI4naC.Contains(testBaseNCBI2naT), false);
    EXPECT_EQ(testBaseNCBI4naG.Contains(testBaseNCBI2naT), false);
    EXPECT_EQ(testBaseNCBI4naR.Contains(testBaseNCBI2naT), false);
    EXPECT_EQ(testBaseNCBI4naS.Contains(testBaseNCBI2naT), false);
    EXPECT_EQ(testBaseNCBI4naM.Contains(testBaseNCBI2naT), false);
    EXPECT_EQ(testBaseNCBI4naV.Contains(testBaseNCBI2naT), false);
}

TEST(AmbiguousBasesTest, TestAmbiguousBaseContainsPureBase)
{
    // Pure Bases
    const char testBaseA{'A'};
    const char testBaseC{'C'};
    const char testBaseG{'G'};
    const char testBaseT{'T'};

    // Diploid Bases
    const char testBaseM{'M'};
    const char testBaseR{'R'};
    const char testBaseS{'S'};
    const char testBaseW{'W'};
    const char testBaseY{'Y'};
    const char testBaseK{'K'};

    // M = A or C
    EXPECT_EQ(ambiguousBaseContainsPureBase(testBaseM, testBaseA), true);
    EXPECT_EQ(ambiguousBaseContainsPureBase(testBaseM, testBaseC), true);
    EXPECT_EQ(ambiguousBaseContainsPureBase(testBaseM, testBaseG), false);
    EXPECT_EQ(ambiguousBaseContainsPureBase(testBaseM, testBaseT), false);

    // R = A or G
    EXPECT_EQ(ambiguousBaseContainsPureBase(testBaseR, testBaseA), true);
    EXPECT_EQ(ambiguousBaseContainsPureBase(testBaseR, testBaseC), false);
    EXPECT_EQ(ambiguousBaseContainsPureBase(testBaseR, testBaseG), true);
    EXPECT_EQ(ambiguousBaseContainsPureBase(testBaseR, testBaseT), false);

    // S = C or G
    EXPECT_EQ(ambiguousBaseContainsPureBase(testBaseS, testBaseA), false);
    EXPECT_EQ(ambiguousBaseContainsPureBase(testBaseS, testBaseC), true);
    EXPECT_EQ(ambiguousBaseContainsPureBase(testBaseS, testBaseG), true);
    EXPECT_EQ(ambiguousBaseContainsPureBase(testBaseS, testBaseT), false);

    // W = A or T
    EXPECT_EQ(ambiguousBaseContainsPureBase(testBaseW, testBaseA), true);
    EXPECT_EQ(ambiguousBaseContainsPureBase(testBaseW, testBaseC), false);
    EXPECT_EQ(ambiguousBaseContainsPureBase(testBaseW, testBaseG), false);
    EXPECT_EQ(ambiguousBaseContainsPureBase(testBaseW, testBaseT), true);

    // Y = C or T
    EXPECT_EQ(ambiguousBaseContainsPureBase(testBaseY, testBaseA), false);
    EXPECT_EQ(ambiguousBaseContainsPureBase(testBaseY, testBaseC), true);
    EXPECT_EQ(ambiguousBaseContainsPureBase(testBaseY, testBaseG), false);
    EXPECT_EQ(ambiguousBaseContainsPureBase(testBaseY, testBaseT), true);

    // K = G or T
    EXPECT_EQ(ambiguousBaseContainsPureBase(testBaseK, testBaseA), false);
    EXPECT_EQ(ambiguousBaseContainsPureBase(testBaseK, testBaseC), false);
    EXPECT_EQ(ambiguousBaseContainsPureBase(testBaseK, testBaseG), true);
    EXPECT_EQ(ambiguousBaseContainsPureBase(testBaseK, testBaseT), true);
}
