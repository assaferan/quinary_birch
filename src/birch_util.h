#ifndef __BIRCH_UTIL_H_
#define __BIRCH_UTIL_H_
// Birch-specific namespace for helper functions

#include "birch.h"

namespace birch_util
{
  template<typename From, typename To>
  To convert_Integer(const From& x);

  int popcnt(Z64 x);

  template<typename From, typename To>
  PrimeSymbol<To> convert_PrimeSymbol(const PrimeSymbol<From>& symbol);

  template<typename From, typename To, size_t n>
  QuadForm<To,n> convert_QuadForm(const QuadForm<From,n>& q);

  template<typename From, typename To, size_t n>
  Isometry<To,n> convert_Isometry(const Isometry<From,n>& s);
  
  template<typename From, typename To, size_t n>
  GenusRep<To,n> convert_GenusRep(const GenusRep<From,n>& from);

  // !! TODO - get rid of this redundant method
  template<typename R>
  R my_pow(const std::map<R,int>& pairs);

  extern int char_vals[256];

  int char_val(W64 x);
}

#include "birch_util.inl"

#endif // __BIRCH_UTIL_H_
