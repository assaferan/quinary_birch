#ifndef FACTORIZATIONS_INCLUDE
#define FACTORIZATIONS_INCLUDE

#include <vector>
#include "gmpxx.h"

std::vector<mpz_class> FactorsInt(mpz_class const& x);
std::vector<mpq_class> FactorsInt(mpq_class const& x);

#endif
