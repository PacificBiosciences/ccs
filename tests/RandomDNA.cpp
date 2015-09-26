
#include "RandomDNA.h"

using namespace std;

string RandomDNA(const size_t n, mt19937* const gen)
{
    constexpr auto bases = "ACGT";
    string result(n, 'A');
    uniform_int_distribution<size_t> rand(0, 3);

    for (size_t i = 0; i < n; ++i)
        result[i] = bases[rand(*gen)];

    return result;
}
