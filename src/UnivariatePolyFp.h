#ifndef __UNIVARIATE_POLY_FP_H_
#define __UNIVARIATE_POLY_FP_H_

#include "birch.h"
#include "Fp.h"
#include "FpElement.h"
#include "UnivariatePoly.h"

template<typename R, typename S>
class UnivariatePolyFp : public UnivariatePoly< FpElement<R,S>, Fp<R,S> >
{
public:
  UnivariatePolyFp(std::shared_ptr< const Fp<R,S>> GF)
    : UnivariatePoly< FpElement<R,S>, Fp<R,S> >(GF) {}

  // create the constant polynomial
  UnivariatePolyFp(const FpElement<R,S> & a) :
    UnivariatePoly< FpElement<R,S>, Fp<R,S> >(a) {}
  
  // create polynomial from coefficients
  UnivariatePolyFp(const std::vector< FpElement<R,S> > & v)
    : UnivariatePoly< FpElement<R,S>, Fp<R,S> >(v) {}

  UnivariatePolyFp(const UnivariatePoly< FpElement<R,S>, Fp<R,S> > & other)
    : UnivariatePoly< FpElement<R,S>, Fp<R,S> >(other) {}

  using UnivariatePoly< FpElement<R,S>, Fp<R,S> >::x;
  
  std::vector< UnivariatePolyFp<R,S> > _sqfFactor(void) const;

  UnivariatePolyInt<R> lift(void) const;

  UnivariatePolyFp<R,S>
  powMod(size_t, const UnivariatePolyFp<R,S> & ) const;

  // assignment and conversion
  UnivariatePolyFp<R,S> & operator=(const UnivariatePoly< FpElement<R,S>, Fp<R,S> > &);
  
  
  // !! TODO -  make it work with inheritance
  static UnivariatePolyFp<R,S> gcd(const UnivariatePolyFp<R,S> & f,
				   const UnivariatePolyFp<R,S> & g);
				    
  static UnivariatePolyFp<R,S> xgcd(const UnivariatePolyFp<R,S> & f,
				    const UnivariatePolyFp<R,S> & g,
				    UnivariatePolyFp<R,S> & s,
				    UnivariatePolyFp<R,S> & t);

  static void divRem(const UnivariatePolyFp<R,S> & f,
		     const UnivariatePolyFp<R,S> & g,
		     UnivariatePolyFp<R,S> & q,
		     UnivariatePolyFp<R,S> & r);
  
  /*
  using UnivariatePoly< FpElement<R,S> >::div_rem;
  using UnivariatePoly< FpElement<R,S> >::gcd;
  using UnivariatePoly< FpElement<R,S> >::xgcd;
  */
  
protected:
  
  std::vector< UnivariatePolyFp<R,S> >
  _czEqDegPartialFactor(size_t r) const;

  std::vector< UnivariatePolyFp<R,S> > _czEqDegFactor(size_t r) const;

  std::vector< UnivariatePolyFp<R,S> > _czDistinctDegFactor() const;
  
};

#include "UnivariatePolyFp.inl"

#endif // __UNIVARIATE_POLY_FP_H_
