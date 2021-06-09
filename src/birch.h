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

/* Builtins */

#define likely(x)   __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

// Constants
constexpr W64 FNV_OFFSET = 0x811c9dc5;
constexpr W64 FNV_PRIME  = 0x01000193;

// Finite fields.
template <typename R, typename S>
class Fp;

typedef Fp<W16,W32>  W16_Fp;
typedef Fp<W32,W64>  W32_Fp;
typedef Fp<W64,W128> W64_Fp;

template <typename R, typename S>
class F2;

typedef F2<W16,W32>  W16_F2;

// Finite field elements
template <typename R, typename S>
class FpElement;

typedef FpElement<W16,W32>  W16_FpElement;
typedef FpElement<W32,W64>  W32_FpElement;
typedef FpElement<W64,W128> W64_FpElement;

// There is no default operator<< for Z128
std::ostream & operator<<(std::ostream & os, const Z128 & z);

#endif // __BIRCH_H_
