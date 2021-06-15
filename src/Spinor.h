#ifndef __SPINOR_H_
#define __SPINOR_H_

#include "birch.h"

template<typename R>
class Spinor
{
  template<typename T>
  friend class Spinor;

public:
  Spinor(const std::vector<R>& primes)
  {
    this->_primes = primes;
    this->_twist = (1LL << this->_primes.size()) - 1;
  }

  template<size_t n>
  inline Z64 norm(const QuadFormZZ<R,n>& q, const Isometry<R,n>& s, const Integer<R>& scalar) const
  {
    Integer<R> tr = Integer<R>::zero();
    for (size_t i = 0; i < n; i++)
      tr += s(i,i);
    // Stub
    // !! TODO - compute the spinor norm of s
    // We should use the genus information for that
    if (n == 3) {
      if (tr != -scalar)
	return this->_computeVals(tr + scalar);
    }
    // for now we let the spinor norm be trivial
    tr = Integer<R>::one();
    return this->_computeVals(tr);
        
  }

  inline const std::vector<R> & primes(void) const
  {
    return this->_primes;
  }

private:
  std::vector<R> _primes;
  Z64 _twist;

  inline Z64 _computeVals(Intger<R> & x) const
  {
    Z64 val = 0;
    Z64 mask = 1;
    for (const R& p : this->_primes)
      {
	while (x % p == 0)
	  {
	    x /= p;
	    val ^= mask;
	  }
	mask <<= 1;
      }
    return val;
  }
};

#endif // __SPINOR_H_
