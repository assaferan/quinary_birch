#include "birch.h"
#include "birch_util.h"
#include "Integer.h"
#include "Rational.h"

namespace birch_util
{
  static int bitcounts[256] = {
			       0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2,
			       3, 3, 4, 3, 4, 4, 5, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3,
			       3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3,
			       4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4,
			       3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5,
			       6, 6, 7, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4,
			       4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5,
			       6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 2, 3, 3, 4, 3, 4, 4, 5,
			       3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 3,
			       4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6,
			       6, 7, 6, 7, 7, 8
  };

  int popcnt(Z64 x)
  {
    if (x <= 0xff) return bitcounts[x];

    int count = 0;
    while (x)
      {
	count += bitcounts[x&0xff];
	x >>= 8;
      }
    return count;
  }

  template<>
  W16 convertInteger<Z>(const Z& x)
  {
    return (W16)mpz_get_ui(x.get_mpz_t());
  }

  template<>
  W16 convertInteger<Z64>(const Z64& x)
  {
    return (W16)x;
  }

  template<>
  W16 convertInteger<Z128>(const Z128& x)
  {
    return (W16)x;
  }
  
  template<>
  W32 convertInteger<Z>(const Z& x)
  {
    return (W32)mpz_get_ui(x.get_mpz_t());
  }

  template<>
  W64 convertInteger<W16>(const W16& x)
  {
    return (W16)x;
  }
  
  template<>
  W32 convertInteger<Z64>(const Z64& x)
  {
    return (W32)x;
  }

  template<>
  W64 convertInteger<Z>(const Z& x)
  {
    return (W64)mpz_get_ui(x.get_mpz_t());
  }

  template<>
  W64 convertInteger<Z64>(const Z64& x)
  {
    return (W64)x;
  }

  template<>
  Z convertInteger<W16>(const W16& x)
  {
    return Z(x);
  }

  template<>
  Z convertInteger<W32>(const W32& x)
  {
    return Z(x);
  }

  template<>
  Z convertInteger<W64>(const W64& x)
  {
    return Z(x);
  }
  
  template<>
  Z convertInteger<Z32>(const Z32& x)
  {
    return Z(x);
  }
  
  template<>
  Z convertInteger<Z64>(const Z64& x)
  {
    return Z(x);
  }

  template<>
  Z convertInteger<W128>(const W128& x)
  {
    W64 top = x >> 64;
    W64 bottom = x & UINT64_MAX;
    
    Z ret = top;
    ret <<= 64;
    ret += bottom;
    return ret;
  }

  template<>
  Z convertInteger<Z128>(const Z128& x)
  {
    W128 abs_x = x < 0 ? -x : x;
    Z ret = convertInteger<W128,Z>(abs_x);
    return x < 0 ? -ret : ret;
  }

  template<>
  Z32 convertInteger<Z>(const Z& x)
  {
    return mpz_get_si(x.get_mpz_t());
  }

  template<>
  Z64 convertInteger<W16>(const W16& x)
  {
    return (Z64)x;
  }

  template<>
  Z64 convertInteger<W32>(const W32& x)
  {
    return (Z64)x;
  }

  template<>
  Z64 convertInteger<W64>(const W64& x)
  {
    return (Z64)x;
  }
  
  template<>
  Z64 convertInteger<Z>(const Z& x)
  {
    return mpz_get_si(x.get_mpz_t());
  }

  template<>
  Z64 convertInteger<Z64>(const Z64& x)
  {
    return x;
  }

  // !! TODO - check that it does what we want
  template<>
  Z128 convertInteger<Z>(const Z& x)
  {
    Z high = x >> 64;
    Z128 res = mpz_get_si(high.get_mpz_t());
    res <<= 64;
    Z low = abs(x) - (high << 64);
    res |= mpz_get_ui(low.get_mpz_t());
    return res;
  }

  template<>
  Z convertInteger<Z>(const Z& x)
  {
    return x;
  }
  
  template<>
  Rational<Z32> convert(const Integer<Z>& x)
  {
    Z32 y = convertInteger<Z,Z32>(x.num());
    return y;
  }

  template<>
  Rational<Z> convert(const Z32 & x)
  {
    Z y = convertInteger<Z32,Z>(x);
    return y;
  }

  template<>
  Rational<Z> convert(const Integer<Z> & x)
  {
    return x;
  }

  template<>
  Rational<Z> convert(const W32 & x)
  {
    Z y = convertInteger<W32,Z>(x);
    return y;
  }
  
  template<>
  Rational<Z> convert(const W64 & x)
  {
    Z y = convertInteger<W64,Z>(x);
    return y;
  }

  int char_vals[256] = {
			1, -1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,
			1, -1,  1, -1, -1,  1,  1, -1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1,
			1, -1, -1,  1,  1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,
			1, -1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1,
			1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,  1, -1, -1,  1,
			1, -1,  1, -1, -1,  1,  1, -1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1,
			1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1,  1, -1, -1,  1, -1,  1,
			1, -1, -1,  1,  1, -1,  1, -1, -1,  1,  1, -1, -1,  1, -1,  1,  1, -1,
			1, -1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1,  1, -1, -1,  1,  1, -1,
			-1,  1, -1,  1,  1, -1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,  1, -1,
			1, -1, -1,  1,  1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,
			1, -1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1,
			1, -1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1,  1, -1, -1,  1,  1, -1,
			-1,  1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1,
			1, -1, -1,  1
  };

}

// There is no default operator<< for Z128
std::ostream & operator<<(std::ostream & os, const Z128 & z)
{
  os << birch_util::convertInteger<Z128,Z>(z);
  return os;
}


