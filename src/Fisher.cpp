// Taken from http://bioinformatics.mpimp-golm.mpg.de/research-projects-publications/supplementary-data/walther/go-term-enrichment-analysis-1/fisher-exact-test-c-code
// and ajdusted to a one-sided test with alternative greater

#include <float.h>
#include <math.h>
#include <pacbio/statistics/Fisher.h>
#include <stdio.h>

namespace PacBio {
namespace Statistics {
double Fisher::fisher_exact_tiss(int chi11, int chi12, int chi21, int chi22)
{
    // int sign = 1;

    int co_occ = chi11;

    const int gene_a = chi11 + chi12;
    const int gene_b = chi11 + chi21;
    const int total_libs = chi11 + chi12 + chi21 + chi22;

    // If the two genes occur few enough times, the minimum number of
    // co-occurrences is 0.  If the total number of times they occur
    // exceeds the number of libraries (say by N), they must overlap
    // at least N times.
    int min_co_occ = 0;
    if (gene_a + gene_b > total_libs) {
        min_co_occ = gene_a + gene_b - total_libs;
    }

    // Maximum number of co-occurrences is at most the number of times
    // the rarer gene occurs in the library :
    int max_co_occ;
    if (gene_a < gene_b) {
        max_co_occ = gene_a;
    } else {
        max_co_occ = gene_b;
    }

    // Calculate the first hypergeometric value

    double base_p = calc_hypergeom(chi11, chi12, chi21, chi22);

    // printf("base_p=%e\n",base_p);

    // If co-occurrences at max possible, then this is our p-value,
    // Also if co-occurrences at min possible, this is our p-value.

    if (co_occ == max_co_occ || co_occ == min_co_occ) {
        // if (co_occ == max_co_occ) {
        //     sign = 1;
        // } else {
        //     sign = -1;
        // }
    } else {
        // Need to add in the other possible p-values.
        double factor_inc = factorInc(chi11, chi12, chi21, chi22);
        double factor_dec = factorDec(chi11, chi12, chi21, chi22);

        // printf("%e,%e:%e\n",factor_inc,factor_dec,base_p);

        // Start out with the current p-value
        double curr_p = base_p;

        // Want to sum the probabilites in the direction of decreasing P
        // if (factor_dec < factor_inc) {
        //     sign = -1;
        //     // Loop down over co-occurrences
        //     do {
        //         // Determine P-value for current chi^2 matrix
        //         curr_p *= factor_dec;

        //         // Add to probability based on recurrence factor
        //         base_p += curr_p;
        //         co_occ--;

        //         // Alter chi^2 matrix to reflect number of co-occurrences
        //         chi11--;
        //         chi22--;
        //         chi12++;
        //         chi21++;

        //         // Get the next value for the recurrence factor
        //         factor_dec = factorDec(chi11, chi12, chi21, chi22);
        //     } while (co_occ > min_co_occ);
        // } else
        if (factor_inc < factor_dec) {
            // sign = 1;
            // Loop up over co-occurrences
            do {
                // Determine P-value for chi^2 matrix from recurrence factor
                curr_p *= factor_inc;

                // Add to probability based on recurrence factor
                base_p += curr_p;
                co_occ++;

                // Alter chi^2 matrix to reflect number of co-occurrences
                chi11++;
                chi22++;
                chi12--;
                chi21--;

                // Get the next value for the recurrence factor
                factor_inc = factorInc(chi11, chi12, chi21, chi22);
                // printf("[%d %d %d %d]\t%e\t%e\t%e\n",chi11, chi12, chi21, chi22, curr_p,base_p,factor_inc);
            } while (co_occ < max_co_occ);
        } else {
            // We are on a saddle point, which means p-value is 1.
            base_p = 1.0;
        }
    }
    return base_p;
}

double Fisher::factorInc(int chi11, int chi12, int chi21, int chi22)
{
    double factor_inc;
    factor_inc = (double)chi12 * chi21;
    factor_inc /= (double)(chi11 + 1) * (chi22 + 1);
    return factor_inc;
}

double Fisher::factorDec(int chi11, int chi12, int chi21, int chi22)
{
    double factor_dec = (double)chi11 * chi22;
    factor_dec /= (double)(chi21 + 1) * (chi12 + 1);
    return factor_dec;
}

double Fisher::gammln(double xx)
{
    static double cof[6] = {76.18009172947146,  -86.50532032941677,    24.01409824083091,
                            -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5};
    double x, tmp, ser;
    int j;

    x = xx - 1.0;

    tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    ser = 1.0;
    for (j = 0; j <= 5; j++) {
        x += 1.0;
        ser += cof[j] / x;
    }
    return -tmp + log(2.50662827465 * ser);
}

double Fisher::factln0(int n)
{
    static double pi = 3.1415926536;

    if (n < 0) {
        //nrerror("Negative factorial in routine FACTLN");
        return 0.0;
    }
    if (n <= 1) return 0.0;
    return 0.5 * log((2.0 * n + 1. / 3.) * pi) + n * (log(1.0 * n) - 1.0);
}

double Fisher::factln(int n)
{
    static double a[101];

    if (n < 0) {
        //nrerror("Negative factorial in routine FACTLN");
        return 0.0;
    }
    if (n <= 1) return 0.0;
    if (n <= 100)
        return a[n] ? a[n] : (a[n] = gammln((double)(n + 1.0)));
    else
        return gammln((double)(n + 1.0));
}

double Fisher::binomialln(int n, int k) { return (factln(n) - factln(k) - factln(n - k)); }

double Fisher::calc_hypergeom(int chi11, int chi12, int chi21, int chi22)
{
    static double b1, b2, b3;
    static int total;

    total = chi11 + chi12 + chi21 + chi22;

    b1 = binomialln(chi11 + chi12, chi11);
    b2 = binomialln(chi21 + chi22, chi21);
    b3 = binomialln(total, chi11 + chi21);

    return exp(b1 + b2 - b3);
}
}
}  //::PacBio::Statistics