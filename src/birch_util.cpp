#include "birch.h"
#include "birch_util.h"

namespace birch_util
{

  template<>
  W16 convert_Integer<Z>(const Z& x)
  {
    return (W16)mpz_get_ui(x.get_mpz_t());
  }

  template<>
  W16 convert_Integer<Z64>(const Z64& x)
  {
    return (W16)x;
  }

  template<>
  W16 convert_Integer<Z128>(const Z128& x)
  {
    return (W16)x;
  }
  
  template<>
  W32 convert_Integer<Z>(const Z& x)
  {
    return (W32)mpz_get_ui(x.get_mpz_t());
  }

  template<>
  W64 convert_Integer<W16>(const W16& x)
  {
    return (W16)x;
  }
  
  template<>
  W32 convert_Integer<Z64>(const Z64& x)
  {
    return (W32)x;
  }

  template<>
  W64 convert_Integer<Z>(const Z& x)
  {
    return (W64)mpz_get_ui(x.get_mpz_t());
  }

  template<>
  W64 convert_Integer<Z64>(const Z64& x)
  {
    return (W64)x;
  }

  template<>
  Z convert_Integer<W16>(const W16& x)
  {
    return Z(x);
  }
  
  template<>
  Z convert_Integer<Z32>(const Z32& x)
  {
    return Z(x);
  }
  
  template<>
  Z convert_Integer<Z64>(const Z64& x)
  {
    return Z(x);
  }

  template<>
  Z convert_Integer<Z128>(const Z128& x)
  {
    W128 abs_x = x < 0 ? -x : x;
    Z ret = convert_Integer<W128, Z>(abs_x);
    return x < 0 ? -ret : ret;
  }

  template<>
  Z convert_Integer<W128>(const W128& x)
  {
    W64 top = x >> 64;
    W64 bottom = x & UINT64_MAX;
    
    Z ret = top;
    ret <<= 64;
    ret += bottom;
    return ret;
  }

  template<>
  Z64 convert_Integer<Z>(const Z& x)
  {
    return mpz_get_si(x.get_mpz_t());
  }

  template<>
  Z64 convert_Integer<Z64>(const Z64& x)
  {
    return x;
  }

  template<>
  Z convert_Integer<Z>(const Z& x)
  {
    return x;
  }

}
