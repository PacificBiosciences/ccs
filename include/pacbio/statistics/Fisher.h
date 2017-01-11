// Taken from http://bioinformatics.mpimp-golm.mpg.de/research-projects-publications/supplementary-data/walther/go-term-enrichment-analysis-1/fisher-exact-test-c-code
/// and ajdusted to a one-sided test with alternative greater
#pragma once

#include <float.h>
#include <math.h>
#include <stdio.h>

namespace PacBio {
namespace Statistics {
class Fisher
{
public:
    static double fisher_exact_tiss(int chi11, int chi12, int chi21, int chi22);

private:
    static double factorInc(int chi11, int chi12, int chi21, int chi22);
    static double factorDec(int chi11, int chi12, int chi21, int chi22);
    static double gammln(double xx);
    static double factln0(int n);
    static double factln(int n);
    static double binomialln(int n, int k);
    static double calc_hypergeom(int chi11, int chi12, int chi21, int chi22);
};
}
}  //::PacBio::Statistics