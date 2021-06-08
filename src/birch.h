#ifndef __BIRCH_H_
#define __BIRCH_H_

#include <gmpxx.h>

/* Type definitions */

typedef mpz_class W;
typedef uint16_t W16;
typedef uint32_t W32;
typedef uint64_t W64;
typedef __uint128_t W128;

typedef mpz_class Z;
typedef int16_t Z16;
typedef int32_t Z32;
typedef int64_t Z64;
typedef __int128_t Z128;

// There is no default operator<< for Z128
std::ostream & operator<<(std::ostream & os, const Z128 & z);

#endif // __BIRCH_H_
