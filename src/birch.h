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

// Where are these used? maybe move them over there?
// Constants
constexpr W64 FNV_OFFSET = 0x811c9dc5;
constexpr W64 FNV_PRIME  = 0x01000193;

// class definitions

// rings (ADT)
template <class Derived, class DerivedElement>
class Ring;

template<class Derived, class DerivedParent>
class RingElement;

// integers
template <typename R>
class IntegerRing;

template <typename R>
class Integer;

// rationals
template<typename R>
class Rational;

template<typename R>
class RationalField;

// finite fields
template<typename R, typename S>
class Fp;

template <typename R, typename S>
class F2;

// Finite field elements
template <typename R, typename S>
class FpElement;

// isometry
template<typename R, size_t n>
class Isometry;

// quadratic forms
template<class R, class Parent, size_t n>
class QuadForm;

template<typename R, typename S, size_t n>
class QuadFormFp;

template<typename R, size_t n>
class QuadFormInt;

template<size_t n>
class Z_QuadForm;

template<size_t n>
class Z64_QuadForm;

template<size_t n>
class Z128_QuadForm;

// matrices
template<class R, class Parent>
class Matrix;

template<class R, class Parent, size_t n>
class SquareMatrix;

// vectors
template<class R, class Parent, size_t n>
class Vector;

// convenient typedefs

// quadratic forms
template<typename R, size_t n>
using R_QuadForm = QuadForm<Integer<R>, IntegerRing<R>, n>;

// matrices
template<typename R, typename S>
using MatrixFp = Matrix< FpElement<R,S>, Fp<R,S> >;

template<typename R>
using MatrixInt = Matrix< Integer<R>, IntegerRing<R> >;

template<typename R>
using MatrixRat = Matrix< Rational<R>, RationalField<R> >;

// square matrices 
template<typename R, typename S, size_t n>
using SquareMatrixFp = SquareMatrix< FpElement<R,S>, Fp<R,S>, n>;

template<typename R, size_t n>
using SquareMatrixInt = SquareMatrix< Integer<R>, IntegerRing<R>, n>;

template<typename R, size_t n>
using SquareMatrixRat = SquareMatrix< Rational<R>, RationalField<R>, n>;

// vectors
template<typename R, typename S, size_t n>
using VectorFp = Vector< FpElement<R,S>, Fp<R,S>, n>;

template<typename R, size_t n>
using VectorInt = Vector< Integer<R>, IntegerRing<R>, n>;

// finite fields
typedef Fp<W16,W32>  W16_Fp;
typedef Fp<W32,W64>  W32_Fp;
typedef Fp<W64,W128> W64_Fp;

typedef F2<W16,W32>  W16_F2;

// finite field elements

typedef FpElement<W16,W32>  W16_FpElement;
typedef FpElement<W32,W64>  W32_FpElement;
typedef FpElement<W64,W128> W64_FpElement;

// Vectors of Finite field elements
template<size_t n>
using W16_VectorFp = VectorFp< W16, W32, n>;
template<size_t n>
using W32_VectorFp = VectorFp< W32, W64, n>;
template<size_t n>
using W64_VectorFp = VectorFp< W64, W128, n>;

// Matrices of finite field elements
typedef MatrixFp< W16, W32> W16_MatrixFp;
typedef MatrixFp< W32, W64> W32_MatrixFp;
typedef MatrixFp< W64, W128> W64_MatrixFp;

/* Struct definitions */

template<typename R>
struct PrimeSymbol {
    R p;
    int power;
    bool ramified;
};

// Prime symbols
typedef PrimeSymbol<Z>   Z_PrimeSymbol;
typedef PrimeSymbol<Z64> Z64_PrimeSymbol;
typedef PrimeSymbol<Z128> Z128_PrimeSymbol;

// There is no default operator<< for Z128
std::ostream & operator<<(std::ostream & os, const Z128 & z);

#endif // __BIRCH_H_
