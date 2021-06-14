#ifndef __UNIVARIATE_POLY_INT_H_
#define __UNIVARIATE_POLY_INT_H_

#include <memory>
#include <set>
#include <unordered_map>
#include <vector>

#include "birch.h"
#include "Integer.h"
#include "UnivariatePoly.h"

template<typename R>
class UnivariatePolyInt : public UnivariatePoly< Integer<R>, IntegerRing<R> >
{
public:
  UnivariatePolyInt()
    : UnivariatePoly< Integer<R>, IntegerRing<R> >(std::make_shared< IntegerRing<R> >())
  {}
  
  UnivariatePolyInt(const Integer<R> & a)
    : UnivariatePoly< Integer<R>, IntegerRing<R> >(a)
  {}

  UnivariatePolyInt(const std::vector< Integer<R> > & v)
    : UnivariatePoly< Integer<R>, IntegerRing<R> >(v)
  {}
  
  UnivariatePolyInt(const UnivariatePoly< Integer<R>, IntegerRing<R> > & other)
    : UnivariatePoly< Integer<R>, IntegerRing<R> >(other)
  {}

  template<typename S, typename T>
  UnivariatePolyFp<S,T> mod(std::shared_ptr< const Fp<S, T> >) const;
  
  std::unordered_map< UnivariatePolyInt<R>, size_t > factor(void) const;

  // hash
  W64 hashValue(void) const;

  // assignment and converstion
  template<class T>
  UnivariatePolyInt<R> & operator=(const UnivariatePolyInt<T> & );
  using UnivariatePoly< Integer<R>, IntegerRing<R> >::operator=;

  // arithmetic
  using UnivariatePoly< Integer<R>, IntegerRing<R> >::operator/;

  // algorithms
  using UnivariatePoly< Integer<R>, IntegerRing<R> >::gcd;

protected:
  // these helper methods are needed for factorization
  
  template<typename S, typename T>
  void _henselStep(std::vector<UnivariatePolyInt<R> > & u,
		   std::vector<UnivariatePolyInt<R> > & v,
		   std::shared_ptr< const Fp<S,T> > GF,
		   size_t i) const;

  template<typename S, typename T>
  std::vector< UnivariatePolyInt<R> >
  _henselLift(const std::vector<UnivariatePolyFp<S,T> > & g,
	      size_t a) const;

  std::vector< UnivariatePolyInt<R> >
  _trialFactor(const std::vector<UnivariatePolyInt<R> > & u,
	       const Integer<R> & N) const;

  void
  _findTrialFactor(const std::vector< UnivariatePolyInt<R> > & u,
		   const Integer<R> & N,
		   size_t & j,
		   std::set<size_t> & C,
		   size_t & s,
		   std::vector< UnivariatePolyInt<R> > & gs);

  Integer<R> _landauMignotte(void) const;

  std::vector< UnivariatePolyInt<R> > _squarefreeFactor(void) const;

  static std::set< std::set<size_t> >
  _subsets(const std::set<size_t> & S, size_t k);
};

namespace std
{
  template<typename R>
  struct hash< UnivariatePolyInt<R> >
  {
    Z64 operator()(const UnivariatePolyInt<R>& p) const
    {
      return p.hashValue();
    }
  };
}

#include "UnivariatePolyInt.inl"

#endif // __UNIVARIATE_POLY_INT_H_
