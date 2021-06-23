#ifndef __BIRCH_UTIL_H_
#define __BIRCH_UTIL_H_
// Birch-specific namespace for helper functions

#include <map>

#include "birch.h"

namespace birch_util
{ 
  
  template<typename From, typename To>
  To convertInteger(const From& x);

  template<typename From, typename To>
  To convert(const From& x);

  int popcnt(Z64 x);

  template<typename From, typename To>
  PrimeSymbol<To> convertPrimeSymbol(const PrimeSymbol<From>& symbol);

  template<typename From, typename To, size_t n>
  QuadFormZZ<To,n> convertQuadForm(const QuadFormZZ<From,n>& q);

  template<typename From, typename To, size_t n>
  Isometry<To,n> convertIsometry(const Isometry<From,n>& s);
  
  template<typename From, typename To, size_t n>
  GenusRep<To,n> convertGenusRep(const GenusRep<From,n>& from);

  template<typename R, size_t n>
  MatrixInt<R> convertMatrix(const SquareMatrixInt<R,n> & from);

  template<typename From, typename To, size_t n>
  SquareMatrixInt<To,n> convertSquareMatrix(const SquareMatrixInt<From,n>& q);
  
  // !! TODO - get rid of this redundant method
  template<typename R>
  R myPow(const std::map<R,int>& pairs);

  template<typename R>
  R gcd(const R &, const R &);

  template<typename R>
  R lcm(const R &, const R &);

  template<typename R>
  R xgcd(const R &, const R &, R &, R&);
  
  int charVal(W64 x);

}

#include "birch_util.inl"

#endif // __BIRCH_UTIL_H_
